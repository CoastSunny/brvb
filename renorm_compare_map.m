function [CM, x] = renorm_compare_map(wh, ico, N, M, varargin)

CM = [];
x = [];
if nargin > 4
    CM = varargin{1};
end

n = N;
m = N+1;
nCM = [];
while n < M
    nfaces = size(ico(n).wh.faces,1);
    mfaces = size(ico(m).wh.faces,1);
    
    if isempty(CM)
%         CM = rand(nfaces/2,3);
%         CM2 = CM;
%         b = [5 6 7 8 1 2 3 4 10 9 12 11 15 16 13 14 17 18 19 20];
%         CM2(b,:) = CM;
%         CM = [CM; CM2];

        CM = rand(nfaces,3);
        b = [25 26 27 28 21 22 23 24 30 29 32 31 35 36 33 34 37 38 39 40];
        CM(b,:) = CM(1:20,:);
%        CM = whitening(CM);

%        b = b+nfaces/2
    end

    f =figure;
    patch('vertices', wh.orig.vertices, 'faces', ico(n).wh.faces, 'FaceVertexCData', linspace(1,nfaces,nfaces)', 'FaceColor', 'flat')
    axis off;
    axis equal;
    cameramenu
%     camlight;
%     lighting phong;
    maxfigsize(f);
    if isempty(nCM)
        colormap(CM);
        nCM = CM;
    else
        colormap(nCM);
    
        f = figure;
        patch('vertices', wh.orig.vertices, 'faces', ico(m).wh.faces, 'FaceVertexCData', linspace(1,mfaces,mfaces)' , 'FaceColor', 'flat') %linspace(1,mfaces,mfaces)'
        axis off;
        axis equal;
        cameramenu
%         camlight;
%         lighting phong;
        colormap(nCM);
        maxfigsize(f);
    end

    f = figure;
    patch('vertices', wh.orig.vertices, 'faces', ico(m).wh.faces, 'FaceVertexCData', ico(m).wh.facesmap , 'FaceColor', 'flat') %linspace(1,mfaces,mfaces)'
    axis off;
    axis equal;
    cameramenu
%     camlight;
%     lighting phong;
    colormap(CM);
    maxfigsize(f);

    nnCM = zeros(nfaces, 3);
    nc = 0;
    vc = 1;
%     nb = zeros(mfaces, 1);
    for i=1:nfaces
        x = find(ico(m).wh.facesmap==i);
%         nnCM(x(1),:) = vc*nCM(i,:) + nc*rand(1,3);
%         nnCM(x(2),:) = vc*nCM(i,:) + nc*rand(1,3);
%         nnCM(x(3),:) = vc*nCM(i,:) + nc*rand(1,3);
%         nnCM(x(4),:) = vc*nCM(i,:) + nc*rand(1,3);

        nnCM(x(1),:) = vc*nCM(i,[3 2 1]);
        nnCM(x(2),:) = vc*nCM(i,[2 1 3]);
        nnCM(x(3),:) = vc*nCM(i,[2 3 1]);
        nnCM(x(4),:) = vc*nCM(i,[3 1 2]);

        mirror = find(b==i);
        if mirror
            y = find(ico(m).wh.facesmap==mirror);
            
            nnCM(y(1),:) = nnCM(x(1),:);
            nb(y(1)) = x(1);
            
            nnCM(y(2),:) = nnCM(x(3),:);
            nb(y(2)) = x(3);
            
            nnCM(y(3),:) = nnCM(x(2),:); 
            nb(y(3)) = x(2);
            
            nnCM(y(4),:) = nnCM(x(4),:);
            nb(y(4)) = x(4);
            
        end
    end
    b = nb;
    
    nCM = nnCM;
%     nCM = whitening(nCM);

    n = n+1;
    m = m+1;
end
end

function CM = whitening(CM)
    MCM = max(CM);
    s = size(CM,1);
    for i=1:s
        CM(i,:) = CM(i,:)./MCM;
    end
end