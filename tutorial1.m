startup_m

%% 0) Initialize
global DATAPATH;
global SUBJ;
global SPMFILES;

SUBJPATH = [DATAPATH filesep SUBJ];
MVBdataFILE = [SUBJPATH filesep 'MVBdata.mat'];
setenv('SUBJPATH', SUBJPATH);

order = 5;

%% 1) Superficies originais - Freesurfer

% cd $FREESURFER_HOME
% source $FREESURFER_HOME/SetUpFreeSurfer.sh
% recon-all -subject subj1 -i $SUBJPATH/dsmri.img -all

% mris_convert subjects/subj1/surf/lh.orig $SUBJPATH/lh.orig.gii
% mris_convert subjects/subj1/surf/rh.orig $SUBJPATH/rh.orig.gii

% mris_convert subjects/subj1/surf/lh.sphere $SUBJPATH/lh.sphere.gii
% mris_convert subjects/subj1/surf/rh.sphere $SUBJPATH/rh.sphere.gii

% mris_convert subjects/subj1/surf/lh.inflated $SUBJPATH/lh.inflated.gii
% mris_convert subjects/subj1/surf/rh.inflated $SUBJPATH/rh.inflated.gii

% mris_convert -c subjects/subj1/surf/lh.curv subjects/subj1/surf/lh.orig $SUBJPATH/lh.curv.gii
% mris_convert -c subjects/subj1/surf/rh.curv subjects/subj1/surf/rh.orig $SUBJPATH/rh.curv.gii

% inside MATLAB

MVBdata.lh.orig = gifti([SUBJPATH filesep 'lh.orig.gii']);
MVBdata.rh.orig = gifti([SUBJPATH filesep 'rh.orig.gii']);

MVBdata.lh.sphere = gifti([SUBJPATH filesep 'lh.sphere.gii']);
MVBdata.rh.sphere = gifti([SUBJPATH filesep 'rh.sphere.gii']);

MVBdata.lh.inflated = gifti([SUBJPATH filesep 'lh.inflated.gii']);
MVBdata.rh.inflated = gifti([SUBJPATH filesep 'rh.inflated.gii']);

MVBdata.lh.curv = gifti([SUBJPATH filesep 'lh.curv.gii']);
MVBdata.rh.curv = gifti([SUBJPATH filesep 'rh.curv.gii']);

save(MVBdataFILE, 'MVBdata');

%% 2) Algoritmo RE - Downsampling

% ATENTIO! Remember to load data into SPM first and create the Forward Model inside 3D source reconstuction!

load(MVBdataFILE, 'MVBdata');

% Use the source localization information from the first file
D = spm_eeg_load([SUBJPATH filesep SPMFILES{1}]);
[MVBdata.ico MVBdata.wh] = renorm_RE(MVBdata.lh, MVBdata.rh, order, D);
%MVBdata.ico = renorm_forward_eeg(MVBdata.ico, order, D);

save(MVBdataFILE, 'MVBdata');

% Now, we are ready to do the source localization and se the results

% 3) and 4) see tutorial2.m
