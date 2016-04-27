#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#define NUM_THREADS 8
#define HALF_PI (1.570796326794897)
#define PI (3.141592653589793)

#include "mex.h"
#include "matrix.h"
#include "lapack_wrapper.h"

/* Input Arguments */

#define FACES_IN            prhs[0]
#define VERTICES_ORIG_IN    prhs[1]
#define VERTICES_SPHERE_IN  prhs[2]
#define FACES_INTER_IN      prhs[3]

/* Output Arguments */

#define FACES_OUT           plhs[0]
#define VERTICES_PER_FACE   plhs[1]
#define AREA_PER_FACE       plhs[2]
#define MEAN_PER_FACE       plhs[3]
#define NORMAL_PER_FACE     plhs[4]
#define FACE_MAP            plhs[5]
#define FACE_MAP_ORI        plhs[6]

#define MEAN_FILE "./debug/mean"
FILE* MEAN;

#define AORI_FILE "./debug/aori"

/* Debug */
#define DEBUG 1

#define HIST_FILE "./debug/hist_vertice_face"
FILE* HIST;

#define INTERMED_FILE_F "./debug/face"
#define INTERMED_FILE_V "./debug/vert"
#define INTERMED_FILE_A "./debug/area"
FILE* INTERMED;

#define STEP_FILE "./debug/step"
FILE* STEP;

/* General Geometric Functions */
typedef double vec3[3];
typedef double vec2[2];

double frand(){return (double)rand()/(double)RAND_MAX;}

double dot(vec3 v1,vec3 v2){ return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];}

void cross(vec3 v1, vec3 v2, vec3* out){
    (*out)[0] = v1[1]*v2[2] - v1[2]*v2[1];
    (*out)[1] = v1[2]*v2[0] - v1[0]*v2[2];
    (*out)[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void diff(vec3 v1, vec3 v2, vec3* out){
    (*out)[0] = v1[0] - v2[0];
    (*out)[1] = v1[1] - v2[1];
    (*out)[2] = v1[2] - v2[2];
}

void sum(vec3 v1, vec3 v2, vec3* out){
    (*out)[0] = v1[0] + v2[0];
    (*out)[1] = v1[1] + v2[1];
    (*out)[2] = v1[2] + v2[2];
}

double length1(vec3 v){  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

double length2(vec3 v1,vec3 v2)
{
    vec3 v;
    v[0] = v2[0] - v1[0]; v[1] = v2[1] - v1[1]; v[2] = v2[2] - v1[2];
    return length1(v);
}

double area(vec3 a, vec3 b, vec3 c){

    vec3 ab, ac, x;

    diff(b,a,&ab);
    diff(c,a,&ac);
    
    cross(ab,ac,&x);

    return .5 * length1(x);

}

/* p * t = a + (b - a) * u + (c - a) * v => if (0 < t,u,v < 1 && (u+v) < 1) ok! */
int pointInTriangle2(vec3 p, vec3 a, vec3 b, vec3 c, int debug){

    double ** A, * cmatrix;
    int i;
    
    A = (double**) calloc(3, sizeof(double*));
    for(i=0; i<3; i++)
        A[i] = (double*) calloc(3, sizeof(double));

    A[0][0] = -p[0];
    A[1][0] = -p[1];
    A[2][0] = -p[2];

    A[0][1] = b[0] - a[0];
    A[1][1] = b[1] - a[1];
    A[2][1] = b[2] - a[2];

    A[0][2] = c[0] - a[0];
    A[1][2] = c[1] - a[1];
    A[2][2] = c[2] - a[2];

    cmatrix = (double *) malloc(9*sizeof(double));

    mat2cvec(3, 3, A, cmatrix);    
    int test2 = matrix_invert(3,cmatrix,cmatrix);
    cvec2mat(3, 3, cmatrix, A); /* Turn it into a proper matrix */
    
    if (!test2){

        double t,u,v;

        t = -(A[0][0] * a[0] + A[0][1] * a[1] + A[0][2] * a[2]);
        u = -(A[1][0] * a[0] + A[1][1] * a[1] + A[1][2] * a[2]);
        v = -(A[2][0] * a[0] + A[2][1] * a[1] + A[2][2] * a[2]);


        test2 = (t >= 0. && t <= 1.) && (u >= 0. && u <= 1.) && (v >= 0. && v <= 1.) && (u+v <= 1.); 

    }else{
        printf ("LAPACK couldn't invert the matrix : %d\n", test2);
        test2 = 0;
    }

    for(i=0;i<3;i++)
        free(A[i]);
    free(A);
    free(cmatrix);

    return test2;
}

struct findpoints_data{
    int beginface;
    int endface;
    int nvertices;
    int * nvout;
    float * aout;
    int ** faces;
    int ** faces_per_vertex;
    int * index_fpv;
    vec3 * verticess;
    vec3 * verticeso;
    float ** polar;
    int ** newfaces;
    int * verticesok;
    int * facesmap_ori;
    int stop;
};


void * findPoints(void * thread_arg){
    int i,j,k, ok=1;
    double l,l2, tmin, tmed, tmax, t[4];
    struct findpoints_data * data;
    data = (struct findpoints_data *) thread_arg; 
    FILE * stop;
    data->stop = 0;

    /*printf("Searching for points from face %d to %d\n", data->beginface, data->endface);*/
    for (i = data->beginface; i < data->endface; i++){
        data->nvout[i] = 0;
        data->aout[i] = 0.0;
        l = 0.0;
        for (j=0; j < data->nvertices; j++){
            
            /* simple test to see if the point is within the minimum and maximum theta of this triangle, and avoid slower precise calculation
            t[4] = data->polar[j][0];
            t[0] = data->polar[data->newfaces[i][0]][0];
            t[1] = data->polar[data->newfaces[i][1]][0];
            t[2] = data->polar[data->newfaces[i][2]][0];

            tmed = t[0]>t[1]?t[0]:t[1];
            tmax = t[2]>tmed?t[2]:tmed;
            tmed = t[0]<t[1]?t[0]:t[1];
            tmin = t[2]<tmed?t[2]:tmed;

            ok = t[4] > tmin && t[4] < tmax;
*/
            if (ok && pointInTriangle2(data->verticess[j], data->verticess[data->newfaces[i][0]], data->verticess[data->newfaces[i][1]], data->verticess[data->newfaces[i][2]], 0)){
                /* simple count */
                data->nvout[i]++; 

                /* area estimation */
                for (k=0; k < data->index_fpv[j]; k++){
                    data->aout[i] += (1.0/3.0) * area( data->verticeso[data->faces[data->faces_per_vertex[j][k]][0]], data->verticeso[data->faces[data->faces_per_vertex[j][k]][1]], data->verticeso[data->faces[data->faces_per_vertex[j][k]][2]]);
                    data->facesmap_ori[data->faces_per_vertex[j][k]] = i;
                }
                
                data->verticesok[j] = 1;
            }
            /* check abort file */
            stop = fopen("./stopsph","r");
            if (stop != NULL) { data->stop = 1; return; }
        }

    }
}

struct findmeans_data{
    int beginface;
    int endface;
    int nvertices;
    int ** faces;
    int ** faces_per_vertex;
    int * index_fpv;
    vec3 * verticess;
    vec3 * verticeso;
    int nnewfaces;
    int ** newfaces;
    float * mean;
    float * normal;
};


void * findMeans(void * thread_arg){
    int i,j,k,l,m,ok;
    double len;
    struct findmeans_data * data;
    data = (struct findmeans_data *) thread_arg; 
    vec3 tmp1, tmp2, norm, mean, normal;

    /*printf("Searching for points from face %d to %d\n", data->beginface, data->endface);*/
    for (i = data->beginface; i < data->endface; i++){
        mean[0] = mean[1] = mean[2] = 0.0;
        normal[0] = normal[1] = normal[2] = 0.0;
        ok=0;
        for (j=0; j < data->nvertices; j++){
            if (pointInTriangle2(data->verticess[j], data->verticess[data->newfaces[i][0]], data->verticess[data->newfaces[i][1]], data->verticess[data->newfaces[i][2]], 0)){

                /* normal estimation */
                for (k=0; k < data->index_fpv[j]; k++){
                        
                    diff(data->verticeso[data->faces[data->faces_per_vertex[j][k]][1]], data->verticeso[data->faces[data->faces_per_vertex[j][k]][0]], &tmp1);
                    diff(data->verticeso[data->faces[data->faces_per_vertex[j][k]][2]], data->verticeso[data->faces[data->faces_per_vertex[j][k]][0]], &tmp2);
                    cross(tmp1, tmp2, &norm);

                    for (m=0; m<3; m++)
                        normal[m] += norm[m];

                }
                
                /* mean estimation */
                for (m=0; m<3; m++)
                    mean[m] += data->verticeso[j][m];
                
                ok++;

            }
        }
        len = length1(normal);
        
        data->mean[i] = mean[0] / (double)ok; 
        data->normal[i] = normal[0] / len;
        
        i+=data->nnewfaces;
        data->mean[i] = mean[1] / (double)ok; 
        data->normal[i] = normal[1] / len;
        
        i+=data->nnewfaces;
        data->mean[i] = mean[2] / (double)ok; 
        data->normal[i] = normal[2] / len;
        
        i-=2*data->nnewfaces;
        
    }
}


double get_coulomb_energy(int N,vec3 p[])
{
    double e = 0;
    int i,j;
    for(i = 0;i<N;i++)  
        for(j = i+1; j<N; j++ ) {
            e += 1/ length2(p[i],p[j]);
        }
    return e;
}

void icosahedron(double R, vec3 * ico, int ** faces){

/*    double phi = (R * R) / sqrt(5.); */
    double phi = (1. + sqrt(5.)) / 2.;

    double t1 = R / sqrt(1. + phi*phi);
    double t = (R*phi) / sqrt(1. + phi*phi);

    ico[0][0]  =  t;    ico[0][1]  =  t1;   ico[0][2]  =  0;
    ico[1][0]  = -t;    ico[1][1]  =  t1;   ico[1][2]  =  0;
    ico[2][0]  =  t;    ico[2][1]  = -t1;   ico[2][2]  =  0;
    ico[3][0]  = -t;    ico[3][1]  = -t1;   ico[3][2]  =  0;

    ico[4][0]  =  t1;   ico[4][1]  =  0; ico[4][2]  =  t;
    ico[5][0]  =  t1;   ico[5][1]  =  0; ico[5][2]  = -t;
    ico[6][0]  = -t1;   ico[6][1]  =  0; ico[6][2]  =  t;
    ico[7][0]  = -t1;   ico[7][1]  =  0; ico[7][2]  = -t;

    ico[8][0]  = 0; ico[8][1]  =  t;    ico[8][2]  =  t1;
    ico[9][0]  = 0; ico[9][1]  = -t;    ico[9][2]  =  t1;
    ico[10][0] = 0; ico[10][1] =  t;    ico[10][2] = -t1;
    ico[11][0] = 0; ico[11][1] = -t;    ico[11][2] = -t1;


/*
    int i,j;
    printf("\nico points\npoints = [");
    for (i=0; i<12; i++){
        for (j=0;j<3;j++)
            printf("%f ", ico[i][j]);
        printf ("; ");

    }
    printf(" ];\n");
        
    for (i=0; i<12; i++)
        printf ("size %d : %f\n", i, length1(ico[i]));
*/

    faces[0][0] =  0; faces[0][1] = 8; faces[0][2] = 4;
    faces[1][0] =  0; faces[1][1] = 5; faces[1][2] = 10;
    faces[2][0] =  2; faces[2][1] = 4; faces[2][2] = 9;
    faces[3][0] =  2; faces[3][1] = 11; faces[3][2] = 5;
    faces[4][0] =  1; faces[4][1] = 6; faces[4][2] = 8;
    faces[5][0] =  1; faces[5][1] = 10; faces[5][2] = 7;
    faces[6][0] =  3; faces[6][1] =  9; faces[6][2] =  6;
    faces[7][0] =  3; faces[7][1] =  7; faces[7][2] =  11;
    faces[8][0] =  0; faces[8][1] =  10; faces[8][2] =  8;
    faces[9][0] =  1; faces[9][1] =  8; faces[9][2] =  10;
    faces[10][0] = 2; faces[10][1] = 9; faces[10][2] = 11;
    faces[11][0] = 3; faces[11][1] = 11; faces[11][2] = 9;
    faces[12][0] = 4; faces[12][1] = 2; faces[12][2] = 0;
    faces[13][0] = 5; faces[13][1] = 0; faces[13][2] = 2;
    faces[14][0] = 6; faces[14][1] = 1; faces[14][2] = 3;
    faces[15][0] = 7; faces[15][1] = 3; faces[15][2] = 1;
    faces[16][0] = 8; faces[16][1] = 6; faces[16][2] = 4;
    faces[17][0] = 9; faces[17][1] = 4; faces[17][2] = 6;
    faces[18][0] = 10; faces[18][1] = 5; faces[18][2] = 7;
    faces[19][0] = 11; faces[19][1] = 7; faces[19][2] = 5;
 
/*
    printf("\nico faces\nfaces = [");
    for (i=0; i<20; i++){
        for (j=0;j<3;j++)
            printf("%d ", faces[i][j]);
        printf ("; ");
    }
    printf(" ];\n");
*/

}

void tetrahedron(double R, vec3 * tetra, int ** faces){

    double a = sqrt(8./3.)*R; 
    
    tetra[0][0] =  a*(sqrt(3.)/3.); tetra[0][1] =  0;       tetra[0][2] = -a*(sqrt(6.)/12.);
    tetra[1][0] = -a*(sqrt(3.)/6.); tetra[1][1] =  a/2.;    tetra[1][2] = -a*(sqrt(6.)/12.);
    tetra[2][0] = -a*(sqrt(3.)/6.); tetra[2][1] = -a/2.;    tetra[2][2] = -a*(sqrt(6.)/12.);
    tetra[3][0] =  0;               tetra[3][1] =  0;       tetra[3][2] =  a*(sqrt(6.)/4.);

    faces[0][0] = 0; faces[0][1] = 1; faces[0][2] = 2;
    faces[1][0] = 0; faces[1][1] = 1; faces[1][2] = 3;
    faces[2][0] = 0; faces[2][1] = 2; faces[2][2] = 3;
    faces[3][0] = 1; faces[3][1] = 2; faces[3][2] = 3;

/* 
    printf("pontos = [%f %f %f;", tetra[0][0],tetra[0][1],tetra[0][2]);
    printf(" %f %f %f;", tetra[1][0], tetra[1][1], tetra[1][2]);
    printf(" %f %f %f;", tetra[2][0], tetra[2][1], tetra[2][2]);
    printf(" %f %f %f]\n", tetra[3][0], tetra[3][1], tetra[3][2]);

    printf("a %f\n", a);
    printf("1 - 2 %f\n", length2(tetra[0],tetra[2]));
    printf("1 - 3 %f\n", length2(tetra[0],tetra[1]));
    printf("1 - 4 %f\n", length2(tetra[0],tetra[3]));
    printf("2 - 3 %f\n", length2(tetra[1],tetra[2]));
    printf("2 - 4 %f\n", length2(tetra[1],tetra[3]));
    printf("3 - 4 %f\n\n", length2(tetra[2],tetra[3]));
*/
}

int findclosest(int nvertices, vec3 * vertices, vec3 newvertice ){

    /* TODO look only in the neighborhood of the vertex */

    /* Find the closest points in the original grid */
    float ll=0.0, llmin=FLT_MAX;
    int closest, k;

    llmin=FLT_MAX;
    for (k=0; k<nvertices; k++){
        ll = length2(newvertice,vertices[k]);
        if (ll < llmin){
            llmin = ll;
            closest = k;
        }
    }
/*  printf("Ponto %d %f %f %f (%f)\n", closest,vertices[closest][0],vertices[closest][1],vertices[closest][2], llmin);  */

    return closest;
}

int getnewindex(int *** edges, int * nextvertice, int * verticemap, int indexori){

    int newindex, i;

    if (!edges[indexori][0][0]){
        edges[indexori][0][0] = indexori;
        edges[indexori][0][1] = (*nextvertice);

        verticemap[(*nextvertice)] = indexori;

        (*nextvertice)++;
    }
    return edges[indexori][0][1];

}

int checkedge(int *** edges, int max, int * nextvertice, int * verticemap, int i, int j, vec3 * verticess, vec3 * tempvertices, float radius){

    /* Sanity Check */
    if (i == j)
        mexErrMsgTxt ("\n\n Face with duplicated Vertex \n\n");

    int master = (i>j)?i:j;
    int slave = (i<j)?i:j;
    int k;

    /* Check to see if this new vertice was already created*/
    for (k = 0; k < max; k++){
        if (!edges[master][k][0]) break;
        if (edges[master][k][0] == slave) return edges[master][k][1];
    }

    if (k < max){ /* Sanity check... */
        /* if not, mark this edge as processed */ 
        edges[master][k][0] = slave;
        edges[master][k][1] = (*nextvertice);

        /* find the midpoint between the two vertices of this edge and update tempvertices */
        tempvertices[(*nextvertice)][0] = verticess[i][0] + .5 * (verticess[j][0] - verticess[i][0]);
        tempvertices[(*nextvertice)][1] = verticess[i][1] + .5 * (verticess[j][1] - verticess[i][1]);
        tempvertices[(*nextvertice)][2] = verticess[i][2] + .5 * (verticess[j][2] - verticess[i][2]);
        double l = length1(tempvertices[(*nextvertice)])/radius;
        tempvertices[(*nextvertice)][0] /= l;
        tempvertices[(*nextvertice)][1] /= l;
        tempvertices[(*nextvertice)][2] /= l;

/*        printf("Distance from new point to 0 : %f to 1: %f (1 to 0 %f) (tah : %f)\n", length2(tempvertices[(*nextvertice)], verticess[i]),length2(tempvertices[(*nextvertice)], verticess[j]),length2(verticess[i], verticess[j]), length1(tempvertices[(*nextvertice)])); */

        verticemap[(*nextvertice)] = -1; /* mark this vertex so I know it is new */
        
        /* update the new vertices count */
        (*nextvertice)++;

        /* return the tempvertices index of this new created vertices */
        return edges[master][k][1];
    }else
        mexErrMsgTxt ("\n\n limit of edges reached \n\n");

}

/* General Convertion and Statistical functions */

void mat2vec(int size, int ** mat, int * vec){
    int l=0, i;
    for (i=0; i<size; i++){
        vec[l] = mat[i][0]; 
        l+=size;
        vec[l] = mat[i][1]; 
        l+=size;
        vec[l] = mat[i][2]; 
        l-=2*size;
        l++;
    }

}

void floatmeanstd(float * x, int xsize, double * mean, double * std){
    double var=0;
    int j;

    (*mean) = 0;
    for (j=0; j < xsize; j++) (*mean) += x[j];
    (*mean) /= (double)xsize;

    for (j=0; j < xsize; j++) var += pow(x[j] - (*mean), 2);
    var /= (double)(xsize-1);
    (*std) = sqrt(var);
 }

void doublemeanstd(double * x, int xsize, double * mean, double * std){
    double var=0;
    int j;

    (*mean) = 0;
    for (j=0; j < xsize; j++) (*mean) += x[j];
    (*mean) /= (double) xsize;

    for (j=0; j < xsize; j++) var += pow(x[j] - (*mean), 2);
    var /= (double) (xsize-1);
    (*std) = sqrt(var);
}

/* Here is the Core of the Program: The Downsampling Algorithm */

void downsample(int finter[], int nverticesinter, int nfacesinter, float voin[], float vsin[], int fin[], int nvertices, int nfaces, int maxfacespervertex, mxArray ** faces_out, mxArray ** vertices_per_face, mxArray ** area_per_face, mxArray ** mean_per_face, mxArray ** normal_per_face, mxArray ** faces_map, mxArray ** faces_map_ori){

    int i, j, k, l, m, n, o, p ;
    float romed, rsmed, romax, rsmax, rsmin, romin, ro, rs;

    time_t inicio, fim;
 


    /* allocate buffers */
    long time1 = clock();
    
    int *counted;
    counted = (int*) calloc(nvertices, sizeof(int));

    int **faces; /* nfaces x 3 */
    faces = (int**) calloc(nfaces, sizeof(int*));
    for (i=0; i<nfaces; i++)
        faces[i] = (int*) calloc(3, sizeof(int));

    vec3 *verticeso; /* nvertices x 3 */
    verticeso = (vec3*) calloc(nvertices, sizeof(vec3));

    vec3 *verticess; /* nvertices x 3 */
    verticess = (vec3*) calloc(nvertices, sizeof(vec3));

    /* polar coordinates for sphere vertexes */
    float **polar; /* nvertices x 2 */
    polar = (float**) calloc(nvertices, sizeof(float*));
    for(i=0; i < nvertices; i++)
        polar[i] = (float*) calloc(2, sizeof(float));

    /* === Copy the original vectors to matrixes === */

    int max=0;
    for (i = 0; i < nfaces; i++){
        /* count the vertices of this face */
        for (j = 0; j < 3; j++){
            faces[i][j] = fin[i]; /* get index of this vertice of this face */
            fin += nfaces;

            if (faces[i][j] >= nvertices){
                printf ("vertice na face %d maior que quantidade de vertices : %d \n", i, faces[i][j]);
                mexErrMsgTxt("There is a face with invalide vertice.");
            }
            
            counted[faces[i][j]]++;
            if (counted[faces[i][j]] > max) 
                max = counted[faces[i][j]];

        }

        fin -= nfaces*3;
    }

    printf ("Maximum number of faces per vertice in the original surface : %d\n", max);

    /* load the vertices in cartesian and polar coordinate */
    double r,x,y;
    for (i = 0; i < nvertices; i++){
        for (j = 0; j < 3; j++){
            verticeso[i][j] = voin[i]; /* get the orininal surface vertice j=[x,y,z] position */
            verticess[i][j] = vsin[i]; /* get the spherical surface vertice j=[x,y,z] position */
            voin += nvertices;
            vsin += nvertices;
        }
        
        r = length1(verticess[i]);
/*        polar[i][0] = acos(verticess[i][2]/r);*/
        polar[i][0] = verticess[i][2]/r;
        x = verticess[i][0];
        y = verticess[i][1];
        polar[i][1] = atan(y/x); 

        if (x > 0)
            polar[i][1] += HALF_PI;
        else
            polar[i][1] += 3. * HALF_PI;

        voin -= nvertices*3;
        vsin -= nvertices*3;
    }

    /* load the faces of each vertice */
    int ** faces_per_vertex;
    faces_per_vertex = (int**) calloc(nvertices, sizeof(int*));
    for (i=0; i<nvertices; i++)
        faces_per_vertex[i] = (int*) calloc(max+1, sizeof(int*));
    
   
    /* and the index for both information */
    int * index_fpv = (int*) calloc(nvertices, sizeof(int));

    double * aori = (double*) calloc(nfaces, sizeof(double*));
    double aoritotal = 0;
    double meanedge = 0;

    for (i = 0; i < nfaces; i++){
        /* count the vertices of this face */
        for (j = 0; j < 3; j++){

            if (faces[i][j] >= nvertices)
                mexErrMsgTxt("There is a face with invalid vertex.");
        
            faces_per_vertex[faces[i][j]][index_fpv[faces[i][j]]] = i;
            index_fpv[faces[i][j]]++;
        }

        /* compute the mean edge*/
        meanedge += length2(verticeso[faces[i][0]], verticeso[faces[i][1]]);
        meanedge += length2(verticeso[faces[i][0]], verticeso[faces[i][2]]);
        meanedge += length2(verticeso[faces[i][1]], verticeso[faces[i][2]]);

        /* compute the area of this face */
        aori[i] = area(verticeso[faces[i][0]], verticeso[faces[i][1]], verticeso[faces[i][2]] );
        aoritotal += aori[i];
    }
/*
    FILE* AORI = fopen(AORI_FILE, "w");
    for (i=0; i < nfaces; i++)
        fprintf(AORI, "%f\n", aori[i]);
    fclose(AORI);
*/

    meanedge /= 3*nfaces;
    printf("Mean vertex distance : %f\n", meanedge);

    /* compute avg and std of original areas */
    double meanoriapf, stdoriapf;
    doublemeanstd(aori, nfaces, &meanoriapf, &stdoriapf);

    printf("Original area : %f [%f (%f)]\n", aoritotal, meanoriapf, stdoriapf );

    /* load the faces from the previous iteration */
    int **facesinter;
    if (finter != NULL){
        facesinter = (int**) calloc(nfaces, sizeof(int*));
        for (i=0; i<nfacesinter; i++)
            facesinter[i] = (int*) calloc(3, sizeof(int));

        for (i = 0; i < nfacesinter; i++){
            for (j=0;j<3;j++){
                facesinter[i][j] = finter[i];
                finter += nfacesinter;
            }
            finter -= nfacesinter*3;
        }
    }

    long time2 = clock();
    printf("#: Time to initialize the buffers  : %lf\n", (time2-time1)/(double)CLOCKS_PER_SEC);
    
    /* === calculate the average radius of the spherical surface === */

    romax = 0.0;
    rsmax = 0.0;
    romin = DBL_MAX;
    rsmin = DBL_MAX;
    romed = 0.0;
    rsmed = 0.0;
    for (i = 0; i < nvertices; i++){
        ro = length1(verticeso[i]);
        rs = length1(verticess[i]);
        if (romax < ro) romax = ro;
        if (rsmax < rs) rsmax = rs;
        if (romin > ro) romin = ro;
        if (rsmin > rs) rsmin = rs;
        romed += ro;
        rsmed += rs;
    }
    /* diferent from romed += ro/nvertices ? sure... */
    romed /= nvertices;
    rsmed /= nvertices;
    printf("Original: Maximum radius %f | Minimum radius %f | Avarage radius %f \n", romax, romin, romed);
    printf("Sphere: Maximum radius %f | Minimum radius %f | Avarage radius %f \n", rsmax, rsmin, rsmed);

    /* === Generate new faces === */

    int nnewvertices;
    int nnewfaces;
    int * facesmap;
    int * facesmap_ori;
    int ** newfaces;
    int ** orinewfaces;
    double meanapf=0, stdapf=0;

    /* temporary new vertices */
    vec3 * tempvertices = (vec3*) calloc (nvertices, sizeof(vec3));

    /* map from new to old vertices */
    int * verticemap = (int*) calloc (nvertices, sizeof(int)); 
    int * oriverticemap = (int*) calloc (nvertices, sizeof(int)); 

    /* map from old to new vertices */
    int * invverticemap = (int*) calloc (nvertices, sizeof(int)); 

    /* faces of each NEW vertice */
    int ** new_faces_per_vertex = (int**) calloc(nvertices, sizeof(int*));
    for (i=0; i<nvertices; i++)
        new_faces_per_vertex[i] = (int*) calloc(max+1, sizeof(int));

    /* and the versors pointing away from the vertice in each face : (B-A) x */
    vec3 ** new_versors_per_vertex = (vec3**) calloc(nvertices, sizeof(vec3*));
    for (i=0; i<nvertices; i++)
        new_versors_per_vertex[i] = (vec3*) calloc(max+1, sizeof(vec3));
 
    int * index_new_fpv = (int*) calloc(nvertices, sizeof(int));

    /* create a face map to original surface */
    (*faces_map_ori) = mxCreateNumericMatrix(nfaces, 1, mxINT32_CLASS, mxREAL);
    facesmap_ori = (int*)mxGetData((*faces_map_ori));

    /* Diferent algorithms for new icosahedron/tetrahedron and for partitioning existing surfaces */
    if (finter == NULL){

        nnewvertices = nverticesinter;

        /* identify the new number of faces */
        if (nnewvertices == 4)
            nnewfaces = 4;
        else if(nnewvertices == 12)
            nnewfaces = 20;
        else
            mexErrMsgTxt("Only Tetrahedron or Icosahedron supported\n");

        /* allocate space for the final set of new faces */
        newfaces = (int**) calloc(nnewfaces, sizeof(int*));
        orinewfaces = (int**) calloc(nnewfaces, sizeof(int*));
        for (i=0; i<nnewfaces; i++){
            newfaces[i] = (int*) calloc(3, sizeof(int));
            orinewfaces[i] = (int*) calloc(3, sizeof(int));
        }

        /* generate "equaly spaced" vertices in the sphere */
        if (nnewfaces == 4)
            tetrahedron(rsmed, tempvertices, newfaces);
        else if(nnewfaces == 20)
            icosahedron(rsmed, tempvertices, newfaces);

        /* Mark all vertices as new so they will be changed by the closest vertices in the spherical surface  */
        for (i=0; i<nnewvertices; i++){
            /* printf ("values for %d : %f %f %f\n", i, tempvertices[i][0],tempvertices[i][1],tempvertices[i][2]); */ 
            verticemap[i] = -1; 
        }

        /* create a dummy face map */
        (*faces_map) = mxCreateNumericMatrix(nnewfaces, 1, mxINT32_CLASS, mxREAL);
        facesmap = (int*)mxGetData((*faces_map));

        for(i=0; i< nnewfaces; i++)
            facesmap[i] = i;

/*
        for (i = 0; i < nnewvertices; i++)
            for (j=i; j< nnewvertices; j++)
                printf("%d - %d %f\n", i, j, length2(verticess[newvertices[i]],verticess[newvertices[j]]));
*/

    }else{

        /* span 4 new faces from each original face */
        
        /* V - E + F = 2 */
        int pnnewvertices = nfacesinter + 2*nverticesinter - 2; /* Eulers Formula - New vertex for each Edge plus old vertices */
        int pnnewfaces = nfacesinter * 4;
        nnewvertices = 0;
        nnewfaces = 0;
    
        printf("pnewvertices %d\n", pnnewvertices);
        printf("pnewfaces %d\n", pnnewfaces);
        printf("nvertices %d\n", nvertices);
        printf("maxfacespervertex %d\n", maxfacespervertex);

        /* final set of new faces */
        newfaces = (int**) calloc(pnnewfaces, sizeof(int*));
        orinewfaces = (int**) calloc(pnnewfaces, sizeof(int*));
        for (i=0; i<pnnewfaces; i++){
            newfaces[i] = (int*) calloc(3, sizeof(int));
            orinewfaces[i] = (int*) calloc(3, sizeof(int));
        }

        /* create the face map (which new face belongs to which old face) */
        (*faces_map) = mxCreateNumericMatrix(pnnewfaces, 1, mxINT32_CLASS, mxREAL);
        facesmap = (int*)mxGetData((*faces_map));

        /* temporary buffers */
        int *** edges = (int***) calloc(nvertices+1, sizeof(int**)); /* map of new vertice (one for each edges) */
        for (i=0; i < nvertices+1; i++){
            edges[i] = (int**) calloc(maxfacespervertex+1, sizeof(int*));
            for (j=0; j < maxfacespervertex+1; j++)
                edges[i][j] = (int*) calloc(2, sizeof(int));
        }
        
        int newpoint01, newpoint02, newpoint12, newindex0, newindex1, newindex2;
        double l;

        for (i=0; i < nfacesinter; i++){

            newindex0 = getnewindex(edges, &nnewvertices, verticemap, facesinter[i][0]);
            newindex1 = getnewindex(edges, &nnewvertices, verticemap, facesinter[i][1]);
            newindex2 = getnewindex(edges, &nnewvertices, verticemap, facesinter[i][2]);

            newpoint01 = checkedge(edges, maxfacespervertex+1, &nnewvertices, verticemap, facesinter[i][0], facesinter[i][1], verticess, tempvertices, rsmed);
            newpoint02 = checkedge(edges, maxfacespervertex+1, &nnewvertices, verticemap, facesinter[i][0], facesinter[i][2], verticess, tempvertices, rsmed);
            newpoint12 = checkedge(edges, maxfacespervertex+1, &nnewvertices, verticemap, facesinter[i][1], facesinter[i][2], verticess, tempvertices, rsmed);

            if (nnewfaces >= pnnewfaces) mexErrMsgTxt("Error predicting number of faces");
            newfaces[nnewfaces][0] = newindex0;
            newfaces[nnewfaces][1] = newpoint01;
            newfaces[nnewfaces][2] = newpoint02;
            facesmap[nnewfaces] = i;
            nnewfaces++;
 
            if (nnewfaces >= pnnewfaces) mexErrMsgTxt("Error predicting number of faces");
            newfaces[nnewfaces][0] = newindex1;
            newfaces[nnewfaces][1] = newpoint01;
            newfaces[nnewfaces][2] = newpoint12;
            facesmap[nnewfaces] = i;
            nnewfaces++;

            if (nnewfaces >= pnnewfaces) mexErrMsgTxt("Error predicting number of faces");
            newfaces[nnewfaces][0] = newindex2;
            newfaces[nnewfaces][1] = newpoint02;
            newfaces[nnewfaces][2] = newpoint12;
            facesmap[nnewfaces] = i;
            nnewfaces++;

            if (nnewfaces >= pnnewfaces) mexErrMsgTxt("Error predicting number of faces");
            newfaces[nnewfaces][0] = newpoint01;
            newfaces[nnewfaces][1] = newpoint02;
            newfaces[nnewfaces][2] = newpoint12;
            facesmap[nnewfaces] = i;
            nnewfaces++;
        }

        /* Free buffers */ 
        for (i=0; i<nvertices; i++){
            for (j=0; j<maxfacespervertex; j++)
                free(edges[i][j]);
            free(edges[i]);
        }
        free(edges);

        printf("pnnewvertices %d => pnnewfaces %d\n", pnnewvertices, pnnewfaces);
    }

    long time3 = clock();
    printf("#: Time to create the new faces : %lf\n", (time3-time2)/(double)CLOCKS_PER_SEC);
    time (&inicio);
    printf ("\n#: BEGIN finding minimal variance distribution : %s\n", asctime (localtime ( &inicio )) );
 
    printf("nnewvertices %d => nnewfaces %d\n", nnewvertices,nnewfaces);

    /* allocate the output buffer */
    (*faces_out) = mxCreateNumericMatrix(nnewfaces, 3, mxINT32_CLASS, mxREAL);
    int * fout = (int*)mxGetData((*faces_out));

    (*vertices_per_face) = mxCreateNumericMatrix(nnewfaces, 1, mxINT32_CLASS, mxREAL);
    int * nvout = (int*)mxGetData((*vertices_per_face));

    (*area_per_face) = mxCreateNumericMatrix(nnewfaces, 1, mxSINGLE_CLASS, mxREAL);
    float * aout = (float*)mxGetData((*area_per_face));

    /* === Now move the newly created vertices until a minimum average number of original vertices per face is reached === */

    /* identify  all vertices that were not counted in any face */
    int * verticesok = calloc(nvertices, sizeof(int));

    /* Make a copy of the original map, so we know which vertices are new */
    for (i=0; i < nnewvertices; i++)
        oriverticemap[i] = verticemap[i];

    /* make a copy of the original new faces to be able to restore original configuration later */
    for (i=0; i< nnewfaces; i++)
        for (j=0; j<3; j++)
            orinewfaces[i][j] = newfaces[i][j];

    pthread_t thread[NUM_THREADS];
    pthread_attr_t attr;
    struct findpoints_data data[NUM_THREADS];

    int step = nnewfaces / NUM_THREADS;
    int begin = 0, rc;

    int freexing = 1;
    n = 0;
    int nmax = 30;
    
    /* initialize the step in direction of the faces with more area then they should have with the mean edge size in the original surface */
    double vstep = 1e-15;
    double minvstep = 100, mindiff = 0.001;
    double laststd = 0.0, bkplaststd = 0.0, lastdiff = 0.0;
    double prevmaxlen = 0.0;
    double aideal = aoritotal / (double) nnewfaces;

    MEAN = fopen(MEAN_FILE, "w");

    /* Allocate memory to save the intermediate state to be able to roll-back in the case of a bad move */
    float * bkp_aout = (float*) mxGetData(mxCreateNumericMatrix(nnewfaces, 1, mxSINGLE_CLASS, mxREAL)); 
    int * bkp_verticemap = (int*) calloc (nvertices, sizeof(int));
    int ** bkp_new_faces_per_vertex = (int**) calloc(nvertices, sizeof(int*));
    for (i=0; i<nvertices; i++)
        bkp_new_faces_per_vertex[i] = (int*) calloc(max+1, sizeof(int));
    vec3 ** bkp_new_versors_per_vertex = (vec3**) calloc(nvertices, sizeof(vec3*));
    for (i=0; i<nvertices; i++)
        bkp_new_versors_per_vertex[i] = (vec3*) calloc(max+1, sizeof(vec3));
    int * bkp_index_new_fpv = (int*) calloc(nvertices, sizeof(int));
    /* NOT FOR THE MOMENT */
    
    do{
        long time4 = clock();

        /* Find the closest points in the cortical sphere and update to it */
        for (i=0; i < nnewvertices; i++){
            if (verticemap[i] < 0)
                verticemap[i] = findclosest(nvertices, verticess, tempvertices[i]); /* then new */
/*                printf ("from %d to %d\n", i, verticemap[i]);*/
        }

        /* update faces */
        for (j=0; j < nnewfaces; j++){
            for (k=0; k<3; k++){
                newfaces[j][k] = verticemap[newfaces[j][k]];
            }
        }

        /* Initialize indexes */
        for (i=0; i<nvertices; i++){
            index_new_fpv[i] = 0;
        }

        /* Map the faces adjacent to the new faces */
        for (i = 0; i < nnewfaces; i++){
            for (j = 0; j < 3; j++){
                new_faces_per_vertex[newfaces[i][j]][index_new_fpv[newfaces[i][j]]] = i;
                index_new_fpv[newfaces[i][j]]++;
            }
        }

        /* For each new vertice, map the versors out from it, pointing in the midle of the face */
        int v1, v2, v3;
        vec3 tmp1, tmp2;
        double len;
        for (i=0; i < nnewvertices; i++)
            for (j = 0; j< index_new_fpv[verticemap[i]]; j++){
                v1 = newfaces[new_faces_per_vertex[verticemap[i]][j]][0];
                v2 = newfaces[new_faces_per_vertex[verticemap[i]][j]][1];
                v3 = newfaces[new_faces_per_vertex[verticemap[i]][j]][2];

                if (v1 == verticemap[i]){
                    diff(verticess[v2], verticess[v1], &tmp1);
                    diff(verticess[v3], verticess[v1], &tmp2);
                }else if (v2 == verticemap[i]) {
                    diff(verticess[v1], verticess[v2], &tmp1);
                    diff(verticess[v3], verticess[v2], &tmp2);
                }else if (v3 == verticemap[i]) {
                    diff(verticess[v1], verticess[v3], &tmp1);
                    diff(verticess[v2], verticess[v3], &tmp2);
                }else mexErrMsgTxt("Algo muito estranho...");

                sum(tmp1, tmp2, &new_versors_per_vertex[verticemap[i]][j]);
                len = length1(new_versors_per_vertex[verticemap[i]][j]);
                new_versors_per_vertex[verticemap[i]][j][0] /= len;
                new_versors_per_vertex[verticemap[i]][j][1] /= len;
                new_versors_per_vertex[verticemap[i]][j][2] /= len;

                len = length1(new_versors_per_vertex[verticemap[i]][j]);

                /*if (i == 1)
                    printf ("(%f) [%d %d] %f %f %f\n", len, j, verticemap[i], new_versors_per_vertex[verticemap[i]][j][0], new_versors_per_vertex[verticemap[i]][j][1], new_versors_per_vertex[verticemap[i]][j][2]);  */

            }

        /* Count how many vertex per face, and the area adjacent to them */

        /* Initialize indexes */
        for (i=0; i<nvertices; i++){
            verticesok[i] = 0;
        }

        /* BEGIN multi-thread */
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        begin = 0;
        for (i = 0; i < NUM_THREADS; i++){
            data[i].beginface = begin;
            if (i+1 == NUM_THREADS)
                data[i].endface = nnewfaces;
            else
                data[i].endface = (begin+step);
            begin = begin + step;        

            data[i].nvertices = nvertices;
            data[i].nvout = nvout;
            data[i].aout = aout;
            data[i].faces = faces;
            data[i].faces_per_vertex = faces_per_vertex;
            data[i].index_fpv = index_fpv;
            data[i].verticess = verticess;
            data[i].verticeso = verticeso;
            data[i].polar = polar;
            data[i].newfaces = newfaces;
            data[i].verticesok = verticesok;
            data[i].facesmap_ori = facesmap_ori;

            rc = pthread_create(&thread[i], &attr, findPoints, (void *)&data[i]); 
            if (rc) mxErrMsgTxt("ERROR; return code from pthread_create() is %d\n", rc);
        }

        pthread_attr_destroy(&attr);
        void * status;
        int stop = 0;
        for(i = 0; i < NUM_THREADS; i++) {
            rc = pthread_join(thread[i], &status);
            if (rc) mxErrMsgTxt("ERROR; return code from pthread_join() is %d\n", rc);
            if (data[i].stop) stop = 1;
        }
        /* END multi-thread */

        /* check to abort file */
        if (stop){
            /* return the last surface */
            break;
        }

        /* write intermediate area in original surface for each new face */
        if (DEBUG){
            char file[255];

            sprintf(file, "%s%d",INTERMED_FILE_F,n);
            INTERMED = fopen(file, "w");
            for (i=0; i< nnewfaces; i++)
                fprintf(INTERMED, "%d %d %d\n", newfaces[i][0]+1,newfaces[i][1]+1,newfaces[i][2]+1);
            fclose(INTERMED);

            sprintf(file, "%s%d",INTERMED_FILE_V,n);
            INTERMED = fopen(file, "w");
            for (i=0; i< nvertices; i++)
                if (oriverticemap[i] != 0 )
                    fprintf(INTERMED, "%f %f %f %d\n", verticeso[verticemap[i]][0],verticeso[verticemap[i]][1],verticeso[verticemap[i]][2], oriverticemap[i]);
            fclose(INTERMED);

            sprintf(file, "%s%d",INTERMED_FILE_A,n);
            INTERMED = fopen(file, "w");
            for (i=0; i< nnewfaces; i++)
                fprintf(INTERMED, "%f\n", aout[i]);
            fclose(INTERMED);
        }

        /* Print all the vertices that werent inside any new face (bug?, in the vertices of the new polihedron?)*/
        printf("\nvertices nao ok\n[");
        int vnok=0;
        for (i=0; i< nvertices; i++)
            if (!verticesok[i]){
                printf("%f %f %f; ", verticess[i][0],verticess[i][1],verticess[i][2]);   
                vnok++;
            }
        printf("]\n %d !!!\n", vnok);

        long time5 = clock();
        printf("\n#: Time to find vertices in triangles : %lf\n", (time5-time4)/(NUM_THREADS*(double)CLOCKS_PER_SEC));

        /* compute average and var */
        floatmeanstd(aout, nnewfaces, &meanapf, &stdapf);

        /* check if the var is decreasing, otherwise reduce step size and go back to previous state 
        if (stdapf - laststd > 0 && laststd) vstep /= 2.0;
        */

        if (stdapf - laststd > 0 && laststd){
            vstep /= 2.0;

            /* RollBack State : verticemap, aout, index_fpv, versors and faces per new vertex
            for (i=0; i < nnewvertices; i++){
                verticemap[i] = bkp_verticemap[i];
                index_new_fpv[verticemap[i]] = bkp_index_new_fpv[verticemap[i]];
                for(k=0; k<index_new_fpv[verticemap[i]]; k++){
                    new_faces_per_vertex[verticemap[i]][k] = bkp_new_faces_per_vertex[verticemap[i]][k];
                    aout[new_faces_per_vertex[verticemap[i]][k]] = bkp_aout[new_faces_per_vertex[verticemap[i]][k]];
                    new_versors_per_vertex[verticemap[i]][k][0] = bkp_new_versors_per_vertex[verticemap[i]][k][0];
                    new_versors_per_vertex[verticemap[i]][k][1] = bkp_new_versors_per_vertex[verticemap[i]][k][1];
                    new_versors_per_vertex[verticemap[i]][k][2] = bkp_new_versors_per_vertex[verticemap[i]][k][2];
                }
            }*/

            printf("Variance increased, rolling back... %f(%f)\n", meanapf, stdapf);
            
            /*
            floatmeanstd(aout, nnewfaces, &meanapf, &stdapf);
            */
        }else{
            /* Save Data 
            for (i=0; i < nnewvertices; i++){
                bkp_verticemap[i] = verticemap[i];
                bkp_index_new_fpv[verticemap[i]] = index_new_fpv[verticemap[i]];
                for(k=0; k<index_new_fpv[verticemap[i]]; k++){
                    bkp_new_faces_per_vertex[verticemap[i]][k] = new_faces_per_vertex[verticemap[i]][k];
                    bkp_aout[new_faces_per_vertex[verticemap[i]][k]] = aout[new_faces_per_vertex[verticemap[i]][k]];
                    bkp_new_versors_per_vertex[verticemap[i]][k][0] = new_versors_per_vertex[verticemap[i]][k][0];
                    bkp_new_versors_per_vertex[verticemap[i]][k][1] = new_versors_per_vertex[verticemap[i]][k][1];
                    bkp_new_versors_per_vertex[verticemap[i]][k][2] = new_versors_per_vertex[verticemap[i]][k][2];
                }
            }
            */
            
            printf("Variance decreased, storing this state...\n");

            n++;
            if (vstep <= minvstep || fabs(laststd-stdapf) <= mindiff)
                freexing = 0;
            else
                freexing = n < nmax;
        }

        printf("Vertices per face #%d : %f (%f) [%f] {%5.15f} <%f> ?%d?\n", n, meanapf, stdapf, aideal, stdapf-laststd, vstep, freexing);
        fprintf(MEAN, "%d %f %f %f %5.15f %f\n", n, meanapf, stdapf, aideal, stdapf-laststd, vstep);

        fflush(stdout);
        mexEvalString("drawnow update;");

        /* Are we going to iterate again? */ 
        if (freexing) {
            
            /* Move the new vertices a step in the direction of the faces with more area then they should have in the original surface */
            double lnmc,lnmc2, dist, dir = 1.0, len = 0.0, maxlen=0.0;

            if (DEBUG){
                char file[255];
                sprintf(file, "%s%d",STEP_FILE,n);
                STEP = fopen(file,"w");
            }
            for (i=0; i < nnewvertices; i++)
                if (oriverticemap[i] < 0){
                    if (index_new_fpv[verticemap[i]]){
                        tempvertices[i][0] = tempvertices[i][1] = tempvertices[i][2] = 0.0; 

                        for(k=0; k<index_new_fpv[verticemap[i]]; k++){
                            len = aout[new_faces_per_vertex[verticemap[i]][k]] - aideal;
/*                            if (fabs(len > maxlen) maxlen = fabs(len);*/
                            if (len > maxlen) maxlen = len;

                            if (prevmaxlen > 0)
                                dir = len / prevmaxlen;
                            else
                                if (len > 0)
                                    dir = 1.0;
                                else
                                    dir = -1.0;

/*                            printf("wtf... %d %f %f %f %f %f\n", k, len, dir*vstep, aout[new_faces_per_vertex[verticemap[i]][k]], aideal); */
                            if (DEBUG)  fprintf(STEP, "%2.10f ", dir*vstep); 

                            tempvertices[i][0] += dir*vstep*new_versors_per_vertex[verticemap[i]][k][0];
                            tempvertices[i][1] += dir*vstep*new_versors_per_vertex[verticemap[i]][k][1];
                            tempvertices[i][2] += dir*vstep*new_versors_per_vertex[verticemap[i]][k][2];
                        }
                        if(DEBUG) fprintf(STEP, "\n"); 
                    }else
                        mxErrMsgTxt("FATAL ERROR! Surrounding faces were not mapped for new vertice %d : %d\n", i, verticemap[i]);

                    tempvertices[i][0] += verticess[verticemap[i]][0]; 
                    tempvertices[i][1] += verticess[verticemap[i]][1]; 
                    tempvertices[i][2] += verticess[verticemap[i]][2]; 

                    lnmc = length1(tempvertices[i]);

                    tempvertices[i][0] /= lnmc;
                    tempvertices[i][1] /= lnmc;
                    tempvertices[i][2] /= lnmc;

                    /* update the vertice back to -1 so it will be replace by the closest point in the cortical sphere in the next iteration */
                    verticemap[i] = -1; 
                }

            if (DEBUG) fclose(STEP);

/*            prevmaxlen = maxlen; */
            prevmaxlen = maxlen/2.0;

            /* Update back the faces to temporary vertices*/
            laststd = stdapf;
            for (j=0; j < nnewfaces; j++)
                for (k=0; k<3; k++)
                    newfaces[j][k] = orinewfaces[j][k];
        }

    }while (freexing);
    
    fclose(MEAN);

    long time6 = clock();
    printf("#: Time to find a best distribution for the new vertices : %lf\n", (time6-time3)/(NUM_THREADS*(double)CLOCKS_PER_SEC));


    time (&fim);
    printf ("\n#: END finding minimal variance distribution : %s\n", asctime (localtime ( &fim )) );

    printf ("\n\nFiding Means and Normals...\n\n");
    fflush(stdout);
    mexEvalString("drawnow update;");

    /* Finish, calculate the mean and the normal of each face */
    (*mean_per_face) = mxCreateNumericMatrix(nnewfaces, 3, mxSINGLE_CLASS, mxREAL);
    float * mean = (float*)mxGetData((*mean_per_face));

    (*normal_per_face) = mxCreateNumericMatrix(nnewfaces, 3, mxSINGLE_CLASS, mxREAL);
    float * normal = (float*)mxGetData((*normal_per_face));
    
    /* BEGIN multi-thread */
    struct findmeans_data dataf[NUM_THREADS];
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    begin = 0;
    for (i = 0; i < NUM_THREADS; i++){
        dataf[i].beginface = begin;
        if (i+1 == NUM_THREADS)
            dataf[i].endface = nnewfaces;
        else
            dataf[i].endface = (begin+step);
        begin = begin + step;        

        dataf[i].nvertices = nvertices;
        dataf[i].faces = faces;
        dataf[i].faces_per_vertex = faces_per_vertex;
        dataf[i].index_fpv = index_fpv;
        dataf[i].verticess = verticess;
        dataf[i].verticeso = verticeso;
        dataf[i].newfaces = newfaces;
        dataf[i].nnewfaces = nnewfaces;
        dataf[i].normal = normal;
        dataf[i].mean = mean;

        rc = pthread_create(&thread[i], &attr, findMeans, (void *)&dataf[i]); 
        if (rc) mxErrMsgTxt("ERROR means; return code from pthread_create() is %d\n", rc);
    }

    pthread_attr_destroy(&attr);
    void * status;
    for(i = 0; i < NUM_THREADS; i++) {
        rc = pthread_join(thread[i], &status);
        if (rc) mxErrMsgTxt("ERROR means; return code from pthread_join() is %d\n", rc);
    }
    /* END multi-thread */

/*
    l = 0;
    double len;
    vec3 tmp1, tmp2, norm;
    for (i=0; i < nnewfaces; i++){
        diff(verticeso[newfaces[i][1]], verticeso[newfaces[i][0]], &tmp1);
        diff(verticeso[newfaces[i][2]], verticeso[newfaces[i][0]], &tmp2);
        cross(tmp1, tmp2, &norm);
        len = length1(norm);

        mean[l] = (1/3.0) * (verticeso[newfaces[i][0]][0] + verticeso[newfaces[i][1]][0] + verticeso[newfaces[i][2]][0]); 
        normal[l] = norm[0] / len;
        l+=nnewfaces;
        mean[l] = (1/3.0) * (verticeso[newfaces[i][0]][1] + verticeso[newfaces[i][1]][1] + verticeso[newfaces[i][2]][1]); 
        normal[l] = norm[1] / len;
        l+=nnewfaces;
        mean[l] = (1/3.0) * (verticeso[newfaces[i][0]][2] + verticeso[newfaces[i][1]][2] + verticeso[newfaces[i][2]][2]); 
        normal[l] = norm[2] / len;
        
        l-=2*nnewfaces;
        l++;
    }
       */ 


    /* output the new faces */
    mat2vec(nnewfaces, newfaces, fout);
 
    /* restore the matlab 1-index based mesh */
    for(i=0; i< nnewfaces*3; i++)
        fout[i]++;
   
    for(i=0; i< nnewfaces; i++)
        facesmap[i]++;

    for(i=0; i< nfaces; i++)
        facesmap_ori[i]++;

    /* === Free buffers === */ 
/*    free(bkp_aout);
    free(bkp_verticemap);
    for (i=0; i<nvertices; i++)
        free(bkp_new_faces_per_vertex[i]);
    free(bkp_new_faces_per_vertex);
    for (i=0; i<nvertices; i++)
        free(bkp_new_versors_per_vertex[i]);
    free(bkp_new_versors_per_vertex);
    free(bkp_index_new_fpv);
*/
    free(orinewfaces);
    free(tempvertices);
    free(verticemap);
    free(oriverticemap);
    free(invverticemap);

    free(counted);

    free(verticesok);

    free(aori);

    for (i=0; i<nfaces; i++)
        free(faces[i]);
    free(faces);

    for (i=0; i<nvertices; i++)
        free(faces_per_vertex[i]);
    free(faces_per_vertex);

    for (i=0; i<nvertices; i++)
        free(new_versors_per_vertex[i]);
    free(new_versors_per_vertex);

    free(index_fpv );

    for (i=0; i<nvertices; i++)
        free(new_faces_per_vertex[i]);
    free(new_faces_per_vertex);

    free(index_new_fpv );
    
    for (i=0; i<nnewfaces; i++)
        free(newfaces[i]);
    free(newfaces);

    free(verticeso);
    
    free(verticess);
    free(polar);

    if (finter != NULL){
        for (i=0; i<nfacesinter; i++)
            free(facesinter[i]);
        free(facesinter);
    }

}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    mwSize nfaces, nvertices, nverticesinter, nfacesinter, t1, t2;

    float *vertices_orig_in_p, *vertices_sphere_in_p;
    int *faces_in_p, *faces_inter_in_p, *faces_out_p, i, j, k, order;

    /* Sanity check for the number of parameters */
    if (nrhs!=4)
        mexErrMsgTxt("Usage: [new_faces, vertices_per_face, area_per_face, means, normals, face_map, face_map_ori] = renorm_invert(faces, sphere_vertices, original_vertices , [ previous_faces | initial_number_of_faces ])");

    /* Sanity check for the dimension of the input vectors */
    nfaces= mxGetM(FACES_IN);
    if (mxGetN(FACES_IN) != 3)
        mexErrMsgTxt("Only Triangle mesh supported for original tesselation");

    printf("Faces %d x 3\n", nfaces);

    nvertices = mxGetM(VERTICES_ORIG_IN);
    if (mxGetN(VERTICES_ORIG_IN) != 3)
        mexErrMsgTxt("Number of dimensions not supported (only 3 dimensions now)");

    printf("Vertices %d x 3\n", nvertices);

    if (mxGetM(VERTICES_SPHERE_IN) != nvertices)
        mexErrMsgTxt("Number of vertices in sphere differs from number of vertices in original vertices of folded surface");

    if (mxGetN(VERTICES_SPHERE_IN) != 3)
        mexErrMsgTxt("Number of dimensions in sphere differs from number of dimentions in original vertices of folded surfaces");

    /* Get the list of vertice in the original surface */
    if (mxIsSingle(VERTICES_ORIG_IN)){
        vertices_orig_in_p = (float*)mxGetData(VERTICES_ORIG_IN);
    }else{
        /* not suported by now 
           checktype(VERTICES_ORIG_IN);
           if (mxIsDouble(VERTICES_ORIG_IN)){
           vertices_orig_in_p = mxGetPr(VERTICES_ORIG_IN);
           }else */
        mxErrMsgTxt("Invalid type for the vertices in the original folded surface (not single)");
    }

    /* Get the list of vertices in the Sphere */
    if (mxIsSingle(VERTICES_SPHERE_IN)){
        vertices_sphere_in_p = (float*)mxGetData(VERTICES_SPHERE_IN);
    }else{
        /* not suported by now 
           checktype(VERTICES_SPHERE_IN);
           if (mxIsDouble(VERTICES_SPHERE_IN)){
           vertices_orig_in_p = mxGetPr(VERTICES_SPHERE_IN);
           }else */
        mxErrMsgTxt("Invalid type for the vertices in the spherical surface (not single)");
    }

    /* Get the list of the original faces */
    faces_in_p = (int*)mxGetData(FACES_IN);
    /* 0-index */
    for (i=0; i<nfaces*3; i++)
        faces_in_p[i]--;

    /* Get the faces in previous iteration (or scalar with initial polyhedra number of faces) */
    nfacesinter = mxGetM(FACES_INTER_IN);

    /* if the 4th parameter ir scalar, try to guess the platonic solid */
    int max = 1;
    if (nfacesinter == 1){

        nfacesinter = floor(mxGetScalar(FACES_INTER_IN));
        if (nfacesinter == 20){
            /* ico  */
            nverticesinter = 12;
            max = 5; 
        }else if(nfacesinter == 4){
            /* tetra */
            nverticesinter = 4;
            max = 3;
        }else
            mexErrMsgTxt("Only Tetrahedron or Icosahedron supported\n");

        faces_inter_in_p = NULL;
    }else{

        if (mxGetN(FACES_INTER_IN) != 3)
            mexErrMsgTxt("Only Triangle mesh supported for intermediary tesselation");

        faces_inter_in_p = (int*)mxGetData(FACES_INTER_IN);

        printf("nfacesinter : %d\n", nfacesinter);

        /* count how many faces per vertex ( and how many vertice ), and correct the faces to 0-index */ 
        nverticesinter = 0;
        int * verticescount = (int*) calloc(nvertices, sizeof(int));
        for (i=0; i<nfacesinter*3; i++){

            /* 0-index */
            faces_inter_in_p[i]--;

            /* update the number of vertices */
            if (!verticescount[faces_inter_in_p[i]])
                nverticesinter++;

            /* update histogram */
            verticescount[faces_inter_in_p[i]]++;
            if (verticescount[faces_inter_in_p[i]] > max)
                max = verticescount[faces_inter_in_p[i]]; 
           
        }

        printf ("nverticesinter : %d\n", nverticesinter);
        printf ("max : %d \n", max);

        if (DEBUG){
            HIST = fopen(HIST_FILE, "w+");
            for (i=0; i<nvertices; i++)
                if (verticescount[i])
                    fprintf(HIST, "%d %d\n", i, verticescount[i]);
            fclose(HIST);
        }

        /* Free buffers */ 
        free(verticescount);
    }

    /* Finally generate the new coarse graned surface */
    time_t inicio, fim;
 
    time (&inicio);
    printf ("Inicio downsample: %s\n", asctime (localtime ( &inicio )) );

    downsample(faces_inter_in_p, nverticesinter, nfacesinter, vertices_orig_in_p, vertices_sphere_in_p, faces_in_p, nvertices, nfaces, max, &FACES_OUT, &VERTICES_PER_FACE, &AREA_PER_FACE, &MEAN_PER_FACE, &NORMAL_PER_FACE, &FACE_MAP, &FACE_MAP_ORI);

    time (&fim);
    printf ("Fim downsample: %s\n", asctime (localtime ( &fim )) );
    
    printf ("\nTotal de tempo : %lf segundos\n", difftime(fim, inicio)); 
    
    /* 1-index */
    if (faces_inter_in_p != NULL)
        for (i=0; i<nfacesinter*3; i++)
            faces_inter_in_p[i]++;

    for (i=0; i<nfaces*3; i++)
        faces_in_p[i]++;

    return;
}
