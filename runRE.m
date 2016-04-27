load(MVBdataFILE, 'MVBdata');

% Use the source localization information from the first file
D = spm_eeg_load([SUBJPATH filesep SPMFILES{1}]);
[MVBdata.ico MVBdata.wh] = renorm_RE(MVBdata.lh, MVBdata.rh, order, D);
%MVBdata.ico = renorm_forward_eeg(MVBdata.ico, order, D);

save(MVBdataFILE, 'MVBdata');


