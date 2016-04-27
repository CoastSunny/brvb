function cfg = renorm_forward(means, normals, varargin)

global SPM8PATH; 

type = 'eeg';
if nargin > 2
    type = varargin{1};
end

% create the set of electrodes
elec = [];
vol = [];

if (nargin > 3)
    D = varargin{2};

%    elec =  D.inv{1}.datareg.sensors;

    scalp = gifti(D.inv{1}.mesh.tess_scalp);
    iskull = gifti(D.inv{1}.mesh.tess_iskull);
    oskull = gifti(D.inv{1}.mesh.tess_oskull);
    sensors = D.inv{1}.datareg.sensors;
    x = sensors.pnt(:,1);
    y = sensors.pnt(:,2);
    z = sensors.pnt(:,3);
    labels = sensors.label;

    vol = ft_read_vol(D.inv{1}.forward.vol);
%     vol.mat = [];
else
    scalp = gifti([SPM8PATH filesep 'canonical' filesep 'scalp_2562.surf.gii']);
    iskull = gifti([SPM8PATH filesep 'canonical' filesep 'iskull_2562.surf.gii']);
    oskull = gifti([SPM8PATH filesep 'canonical' filesep 'oskull_2562.surf.gii']);

    [labels x y z] = textread([SPM8PATH filesep 'EEGtemplates' filesep 'biosemi128.sfp'], '%s %f %f %f');
%    [labels x y z] = textread([SUBJPATH filesep 'elec.pos'], '%s %f %f %f');
    
    % align
    x = x - double(mean(x)) + double(mean(scalp.vertices(:,1)));
    y = y - double(mean(y)) + double(mean(scalp.vertices(:,2)));
    z = z - double(max(z)) + double(max(scalp.vertices(:,3)));
end

% remove reference electrodes
elec.pnt = [x y z];
elec.label = labels;
i=1;
while (i <= size(elec.label,1))
    if ~isempty(strfind(elec.label{i},'nas')) || ~isempty(strfind(elec.label{i},'lpa')) || ~isempty(strfind(elec.label{i},'rpa'))
        elec.pnt(i,:) = [];
        elec.label(i,:) = [];
    else
        i = i+1;
    end
end

if (~isempty(vol))

    fprintf('Calculating the volume conduction using BEM model from SPM\n');
    [cfg.vol, cfg.elec] = ft_prepare_vol_sens(vol, elec, 'channel', D.inv{1}.forward.channels);
    
else
    if strcmp(type , 'eeg')

        % create volume conduction model for EEG (3 concentric spheres)
        headshape(1).pnt = double(iskull.vertices);
        headshape(2).pnt = double(oskull.vertices);
        headshape(3).pnt = double(scalp.vertices);

        cfg.headshape    = headshape;
        cfg.conductivity = [1 1/80 1];

        [vol, cfg] = ft_prepare_concentricspheres(cfg);

        % project the electrodes in the outmost sphere
        [cfg.vol, cfg.elec] = ft_prepare_vol_sens(vol, elec);

    else
        if strcmp(type , 'meg')

            % create the volume conduction model for MEG (1 sphere)
            cfg.headshape    = double(scalp.vertices);
            [vol, cfg] = ft_prepare_singleshell(cfg);

            elec.type = 'meg';
            elec.pnt = elec.pnt * 1.2;  % scale them to a unit sphere and shift outward a bit
            elec.ori = elec.pnt;        % unit length
            ne = size(elec.pnt,1);
            for i=1:ne
              elec.label{i} = sprintf('meg%03d', i);
            end
            elec.tra = eye(ne, ne);

            % project the electrodes in the outmost sphere
            [cfg.vol, cfg.elec] = ft_prepare_vol_sens(vol, elec);
        else
            vol.mat = eye(size(elec.pnt,1));
            cfg.vol = vol;
            cfg.elec = elec;
        end
    end
end

ndips = size(normals,1);
% different versions of fieldtrip
if isfield(cfg.elec, 'pnt')
    nsens = size(cfg.elec.pnt,1);
else
    nsens = size(cfg.elec.chanpos,1);
end

% Start parallel processing, if necessary
isOpen = matlabpool('size') > 0;
if ~isOpen
   fprintf('Oppening matlabpool...\n');
   matlabpool;
end
        
% compute the leadfield for a dipole at each face (extended)
cfg.dip.pos = means;
cfg.dip.mom = normals;
cfg.lf = ft_compute_leadfield(cfg.dip.pos, cfg.elec, cfg.vol);

% compute the leadfield along the normals (reduced or dipole fixed orientation)
cfg.rlf = zeros(nsens, ndips);
for i = 1:ndips
    cfg.rlf(:, i) = cfg.lf(:, (3*i - 2):(3*i)) * normals(i, :)';
end
