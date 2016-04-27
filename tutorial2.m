
%% 0) Initialize

startup_m
global DATAPATH;
global SUBJ;
global SPMFILES;

SUBJPATH = [DATAPATH filesep SUBJ];
load([SUBJPATH filesep 'MVBdata.mat'], 'MVBdata');

order = 5;

%% 3) Source Localization

curpath = pwd;
cd(SUBJPATH);

% V = GJ
MaxVBIterations = 2000;
MaxEps = 1.0e-5;

% MVBdata.ico(order).cfg = renorm_simulate3(MVBdata.wh.orig.vertices, MVBdata.ico(order).wh.faces, MVBdata.ico(order), false, 15, 10, true);

nfiles = length(SPMFILES);
states = cell(1,nfiles);
k = 1;
for j = 1:nfiles
    D = spm_eeg_load([SUBJPATH filesep SPMFILES{j}]);
    for i = 1:D.ntrials
        
        % load removing pre-stimulus time
        V = D(D.meegchannels,find(D.time>=0),i);
        % scale from micro-volts to volts
        V = V*1e-6;
        % convert to double to use Sparse matrixes
        V = double(V);

        % initialize CFG with the data
        cfg = MVBdata.ico(order).cfg;
        % the potential
        cfg.rpot = V;
        % and real dipole (for real data its zeros (used for simulation))
        cfg.dip.rmom = zeros(size(cfg.dip.mom,1), size(V,2));

        % save this configuration for this state/trial
        states{k}.cfg = cfg;
        states{k}.label = D.conditions{i};
        k = k + 1;
    end
end

% run the actual source localization for each trial/state
renorm_batch(MVBdata.wh, MVBdata.ico, 1, order, 0, size(states,2), MaxVBIterations, MaxEps, false, states, true);

cd(curpath);

%% 4) Visualização

% create real space with the same amount of time for plotting function
result = renorm_load([SUBJPATH filesep 'batch' filesep], 100, true);

initialimiar = 100;
invertealpha = true;
alphas = false;

renorm_plot(MVBdata.wh.orig, MVBdata.ico, result{1}, initialimiar, invertealpha, order, alphas);
% renorm_plot_ori(MVBdata.wh, MVBdata.ico, result{1}, initialimiar, invertealpha, order, alphas);
