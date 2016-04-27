function renorm_batch_ori(wh, ico, mo, o, NSR, varargin)

global DATAPATH;
global SUBJ;
global ANALISES;

batchdir = [DATAPATH filesep SUBJ filesep ANALISES ];
try
    mkdir(batchdir);
catch
end

global DEBUG

files = dir(batchdir);
inin = 0;
for i=1:size(files,1)
    if str2num(files(i).name) %#ok<ST2NM>
        inin = inin+1;
    end
end

N = 10;
if (nargin > 5)
    N = varargin{1};
end
R = 10;
if (nargin > 6)
    R = varargin{2};
end
meps = 1.0e-6;
if (nargin > 7)
    meps = varargin{3};
end
sim=1;
if (nargin > 8)
    sim = varargin{4};
end
limsup=10;
liminf=5;
smear=true;
oldstates = [];
orisig = false;
if (nargin > 9)
    if sim
        limsup=varargin{5};
        if (nargin > 10)
            liminf=varargin{6};
        end
        if (nargin > 11)
            smear=varargin{7};
        end
    else
        oldstates = varargin{5};
        if (nargin > 10)
            orisig=varargin{6};
%             if (nargin > 11)
%                 SGE = varargin{7};
%             end
        end
%         fprintf('nargin = %d %d\ n', nargin, orisig);
    end
end


%RB = cell(1,N); %preallocate the cell array
%RB{1} = struct('dip',[],'dipori',[], 'dipmmq', [], 'SS', [], 'SSOri', [], 'cfg', []);

STATE = struct('dip',[],'dipori',[], 'dipmmq', [], 'SS', [], 'SSOri', [], 'cfg', [], 'X', []);

ini=tic;
for n=1:N
    t = tic;
    fprintf('\n\nTrial %d\n\n', n+inin);

    % run another simulation
    cfg=[];
    label = num2str(n+inin);
%     SGE = 1;
    if sim
        if sim == 1
            fcfg = renorm_simulate3(wh.orig.vertices, ico(o).wh.faces, ico(o), false, limsup, liminf, smear);
        elseif sim == 2
            fcfg = renorm_simulate_orig(wh.orig.vertices, wh.orig.faces, wh.cfg, false, limsup, liminf);
        else
            error('Unknow simulation type.');
        end
        
        % keep only relevant information (otherzise cfg is too big to store)
        X = fcfg.rpot; % corrupted field (will be updated)
        cfg.dip.hot = fcfg.dip.hot; % orignal current
        cfg.dip.rmom = fcfg.dip.rmom; % orignal current (distributed)
        cfg.rpot = fcfg.rpot; % original field
    else
        if (~isempty(oldstates))
            fprintf('OK! aproveitando simulação %d ... \n\n', n);
            pause(1);

            if isfield(oldstates{n}, 'STATE')
                state = oldstates{n}.STATE;
            else
                state = oldstates{n};
            end
            
            cfg = state.cfg;
            if isfield(state, 'label')
                label = strcat(state.label, '-', label);
            end
            if orisig
                X = cfg.rpot;
                fprintf('Campo original, corrompendo com ruido NSR = %f .\n', NSR);
                pause (1);
            else
                X = state.X;
                NSR = 0;
                fprintf('Campo ja corrompido, removendo ruido.\n');
            end
            
        else
            fprintf('Warning! Not simulating new scenario every trial... \n\n');
        end
    end

    STATE.cfg = cfg;
%    [STATE.dip, STATE.dipori, STATE.dipmmq, STATE.SS, STATE.SSOri] = renorm_invert3(ico, cfg.rpot, R);
%     figure; plot(X(62,:), '.r');
%     drawnow
    [STATE.dip, STATE.dipori, STATE.dipmmq, STATE.SS, STATE.SSOri, STATE.X] = renorm_invert13(ico, mo, o, X, R, meps, NSR);

    
    thisbatchdir = strcat(batchdir, filesep, label);
    mkdir (thisbatchdir);
    if DEBUG
        movefile([batchdir filesep 'S*.dat'], thisbatchdir);
        movefile([batchdir filesep 'J*.dat'], thisbatchdir);
    end
    save([thisbatchdir filesep 'STATE.mat'], 'STATE');

    ttrial = toc(t);
    fprintf('\nTrail %d : %f \n', n, ttrial );
end

ttotal = toc(ini);
fprintf ('\n\nTotal Batch Time : %f \n', ttotal);
