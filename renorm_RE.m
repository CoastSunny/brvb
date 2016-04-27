function [ico wh] = renorm_RE(lh, rh, o, varargin)

global DATAPATH;

debugdir = [pwd filesep 'debug' filesep];
try
    mkdir(debugdir);
catch
end

forward = 1; % default EEG spherical
if (nargin > 3) % 2 is MEG spherical or the SPM D variavel
    forward = varargin{1};
end

printcompare = false;
if (nargin > 4) 
    printcompare = true;
end


% Iterate RE Algorithm
icotmpfile = fullfile(DATAPATH,'icotmp.mat');
if exist(icotmpfile, 'file')
    x = icotmpfile;
    ico = x.ico;
else
    ico = [];
end

i=length(ico);

while i < o

    i = i+1;

    % NEW, GENRATE MAP ORI
    if i == 1
        % first iretation of COSURE must be the ICOSAHEDRON (or thetahedron)
        [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap, ico(i).rh.facesmapori] = renorm_downsample_new(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, 20);
%         [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap, ico(i).rh.facesmapori] = renorm_downsample_new_unbalanced(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, 20);
        saveintermed('rh', debugdir, i);
        [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap, ico(i).lh.facesmapori] = renorm_downsample_new(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, 20);
%         [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap, ico(i).lh.facesmapori] = renorm_downsample_new_unbalanced(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, 20);
        saveintermed('lh', debugdir, i);
    else
        % second on is the biossection of the previous surface
        [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap, ico(i).rh.facesmapori] = renorm_downsample_new(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, ico(i-1).rh.faces);
%         [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap, ico(i).rh.facesmapori] = renorm_downsample_new_unbalanced(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, ico(i-1).rh.faces);
        saveintermed('rh', debugdir, i);
        [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap, ico(i).lh.facesmapori] = renorm_downsample_new(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, ico(i-1).lh.faces);
%         [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap, ico(i).lh.facesmapori] = renorm_downsample_new_unbalanced(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, ico(i-1).lh.faces);
        saveintermed('lh', debugdir, i);
    end
    
    % ORIGINAL, NO MAP ORI
%     if i == 1
%         % first iretation of COSURE must be the ICOSAHEDRON (or thetahedron)
%         [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap] = renorm_downsample(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, 20);
%         saveintermed('rh', debugdir, i);
%         [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap] = renorm_downsample(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, 20);
%         saveintermed('lh', debugdir, i);
%     else
%         % second on is the biossection of the previous surface
%         [ico(i).rh.faces, ico(i).rh.vperface, ico(i).rh.aperface, ico(i).rh.mean, ico(i).rh.normal, ico(i).rh.facesmap] = renorm_downsample(rh.orig.faces, rh.orig.vertices, rh.sphere.vertices, ico(i-1).rh.faces);
%         saveintermed('rh', debugdir, i);
%         [ico(i).lh.faces, ico(i).lh.vperface, ico(i).lh.aperface, ico(i).lh.mean, ico(i).lh.normal, ico(i).lh.facesmap] = renorm_downsample(lh.orig.faces, lh.orig.vertices, lh.sphere.vertices, ico(i-1).lh.faces);
%         saveintermed('lh', debugdir, i);
%     end
    
    save(icotmpfile, 'ico');

end

% join both hemispheres
[ico wh] = renorm_join(ico, rh, lh);

% calculate average distance between vertices
ico = renorm_meandist(ico, o);

% calculate average face areas
wh.areas = estimate_area(wh,ico);    

% find nearst neighbors
for i=1:o
    [ico(i).wh.vizinhos ico(i).wh.vizinhos2]= renorm_vizinhos(ico(i).wh.faces);
end

% create the compare map (just to visualize results)
if printcompare
    for i=1:o-1
        renorm_compare_map(wh, ico, i, i+1);

    end
end

% compute the Lead Field for each surface
if isnumeric(forward)
    if forward == 1
        ico = renorm_forward_eeg(ico, o);
    else
        if forward == 2
            ico = renorm_forward_meg(ico, o);
        end
    end
else
    ico = renorm_forward_eeg(ico, o, forward);
end

end

function saveintermed(h, debugdir, i)
    dir = [debugdir 'ico' h num2str(i) filesep];
    if ~exist(dir, 'dir')
        mkdir (dir);
    end
        movefile([debugdir 'face*'], dir, 'f');
        movefile([debugdir 'area*'], dir, 'f');
        movefile([debugdir 'vert*'], dir, 'f');
        movefile([debugdir 'mean'], dir, 'f');
    try
        movefile([debugdir 'step*'], dir, 'f');
    catch exep
        fprintf('No movements to balance vertex necessary (unbalanced grids?) %s\n',exep.message);
    end
end


%     else
%         if forward == 3
%             printf('Calculating original LeadField for ico = %d\n', o);
%             ico(o).cfg = renorm_forward(ico(o).wh.mean, ico(o).wh.normal);
%             i = o-1;
%             while i >= 1
%                 ico(i).cfg = ico(i+1).cfg;
%                 
%                 % compute the leadfield as an average of the next order
%                 % leadfield
%                 ico(i).cfg.dip.pos =ico(i).wh.mean;
%                 ico(i).cfg.dip.mom = ico(i).wh.normal;
%                 ndips = size(ico(i).wh.mean,1);
%                 nsens = size(ico(i).cfg.elec.pnt,1);
%                 
%                 ico(i).cfg.lf = zeros(nsens,ndips*3);
%                 ico(i).cfg.rlf = zeros(nsens,ndips);
%                 
%                 printf('Calculating average LeadField for ico = %d\n', i);
%                 size(ico(i).cfg.rlf )
%                 for j=1:ndips
%                     thisj = find(ico(i+1).wh.facesmap == j);
%                     ico(i).cfg.lf(:,(3*j)-2) = sum(ico(i+1).cfg.lf(:,(3*thisj)-2),2)./size(thisj,1);
%                     ico(i).cfg.lf(:,(3*j)-1) = sum(ico(i+1).cfg.lf(:,(3*thisj)-1),2)./size(thisj,1);
%                     ico(i).cfg.lf(:,(3*j)) = sum(ico(i+1).cfg.lf(:,3*thisj),2)./size(thisj,1);
%                     
%                     ico(i).cfg.rlf(:, j) = sum(ico(i+1).cfg.rlf(:,thisj),2)./size(thisj,1);
%                 end
%                 size(ico(i).cfg.rlf )
%                 
%                 i = i-1;
%                 
%             end
%         end
