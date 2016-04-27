function [J, Jori, Jmmq, SS, Sori, X] = renorm_invert13(ico, mino, finer, X, n, meps, varargin)

% This function invert the solution of the EEG problem
% It solves the "inverse problem" maximizing the Free Energy (Sato et al
% 2004)
% 
% Different from the original algorithm, this function uses the information
%   of coarse grained versions of the cortex (ico(i) , i = 1:n) as prior
%   information in the finner versions
%
% Parameters:
%
% ico : the icosahedron aproximations and the original source space
% ico(i).wh.facesmap : maps the faces in this scale to the faces in the
% previous scale
% ico(i).cfg.rlf : the leadfielde in this scale
%
% X : MxT matrix, where M is the number of sensors (rows) and T is the
% time, or samples of each sensor (columns)
%
% n : maximum number of VB iterations
%

global DATAPATH;
global SUBJ;
global ANALISES;

batchdir = [DATAPATH filesep SUBJ filesep ANALISES ];

global DEBUG;

% Number of Channels
Nsens = size(X,1);

% Duration of the signal
T = size(X, 2);
HT = ceil(T/2);

% Add gaussian noise
NSR=0;
if (nargin > 6)
    NSR = varargin{1};
end
    
if NSR
%        X2 = X'*X/Nsens
%        X2 = mean(X,2);
%        X2 = mean(abs(X),2);
    Xm = sqrt(X'*X/Nsens);
    for i=1:size(Xm,2)
        X(:, i) = X(:,i) + sqrt(NSR) * Xm(i,i) * randn(size(X,1), 1);
    end
end

fprintf ('NSR = %f\n', NSR);

% Skipping first few coarse grained surfaces?
%mino = 2;

% Finer scale, the one that "generated the measured field"
% finer = size(ico,2);
% while (finer > 0 && (~isfield(ico(finer), 'cfg') || ~isfield(ico(finer).cfg, 'rlf')))
%     finer = finer-1;
% end
% 
% if (finer == 0)
%     warning('Renorm:nofinerscale', 'No scale has forward solution');
%     return
% end
%finer = 5;

fprintf('finer = %d\n\n', finer)

% Before everything else, set initial values and convergence criterea
%alphalim = 80;
histcut = 2;
treshold = 100;
maxsnfaces = n;
scale=finer-mino;

a0 = 1;
scalea = 1;
ta0 = a0 * scalea^scale;

g0 = 0.1;
%g0 = 1;
scaleg = 1;
tg0 = g0 * scaleg^scale;

% Noise ariance per Channel
% Estimate noise covariance
%error('a')
%SGini = 1;
%SG = 1e22*speye(Nsens, Nsens);
%SG = NSR*1e5*speye(Nsens, Nsens);
%SG = speye(Nsens, Nsens);
% mean(diag(SG))

XT= sqrt(sum(X.^2,2));
%SGE = sqrt(1/(mean(XT.^2) - mean(XT)^2))
SGE = 1/(mean(XT.^2) - mean(XT)^2)
SG =  SGE * speye(Nsens, Nsens);

%% First, run the algorithm to the finer scale only to compare with multi-scale
ostart=tic;

G = ico(finer).cfg.rlf;
NJs = size(G,2);

Aori = a0*speye(NJs, NJs);
A0ori = Aori;
Jori=[];

% % estimate the noise
%[~, ~, SGini] = renorm_VB_sparse8(X, G, Aori, A0ori, g0, SG);

if DEBUG
    fAid = fopen([batchdir filesep 'Sori.dat'], 'w');
    fJid = fopen([batchdir filesep 'Jori.dat'], 'w');
end

prevmalpha = 0;
neg = false;
snfaces = 0;
prevnfaces = 0;
prevAcond = 0;
betaori = 0;
stopfile = 0;

for i=1:n
    % test stop file
    stopfile = fopen('stopinv');
    if (stopfile > 0)
        warning('Stoping because of stopfile!'); %#ok<WNTAG>
        break;
    end
    
    % write debug
    if DEBUG
        fprintf(fAid, '%3.10f ', full(diag(Aori)));
        fprintf(fAid, '\n');
    end
    
    % save working state
    AoriP = Aori;
    JoriP = Jori;
    betaoriP = betaori;
 
    % Try new iteration
    [Jori, Aori, betaori] = renorm_VB_sparse8(X, G, Aori, A0ori, g0, SG);

    % write debug
    if DEBUG
        fprintf(fJid, '%3.10f ', full(Jori(:,HT)));
        fprintf(fJid, '\n');
    end

    DA = full(diag(Aori));
%    malpha = mean(DA);
    ga = sort(DA);
    malpha = mean(ga(1:ceil(NJs*0.1)));
    maxalpha = max(DA);
    minalpha = min(DA);
    Acond = abs(log10(condest(Aori)));
    Acond1mom = abs(Acond - prevAcond);
    fprintf('VB Only step %d (%1.1e %1.1e)\t', i, a0, g0);
    fprintf('Alpha : > %e %e %e < | ', malpha, maxalpha, minalpha);
    fprintf('Beta : %e %e \n', betaori, abs(betaori-betaoriP));

    % Check for numerical errors
    if (minalpha < 0)
        Aori = AoriP;
        Jori = JoriP;
        neg = 1;
        fprintf('WARNING! alpha assumed negative values! Keeping previous calculation\n');
        break;
    end
    
    % Check for stop criterea
    aeps = abs(prevmalpha-malpha);
    eps = aeps/abs(malpha);
    [jm, ~] = hist(abs(mean(Jori,2)), treshold);
    nfaces = sum(jm(histcut:end),2);
    if (nfaces == prevnfaces)
        snfaces = snfaces+1;
    else
        snfaces = 0;
    end
    
    maxJ = full(mean(max(abs(Jori))));
    meanJ = full(mean(mean(abs(Jori))));
           
    fprintf ('Alpha change: %e (%1.0e) ', eps, meps);
    fprintf ('| J (%f) : %d (%d) {%d} | [%f (%f)] ', 1-histcut/treshold, nfaces, snfaces, ceil(maxsnfaces), maxJ, meanJ);
    fprintf ('| Acond : {%f (15) (%e)} \n', Acond, Acond1mom);

    if (prevmalpha > 0 && i > 100 && eps < meps) % || Acond > 15)
        break;
    end
    
    prevmalpha = malpha;
    prevnfaces = nfaces;
    prevAcond = Acond;
end
if (i == n)
    fprintf('Original iteration exit for maximum iteration. Considering change minimum alpha variation\n');
end

if DEBUG
    fclose(fAid);
    fclose(fJid);
end

oend = toc(ostart);
fprintf('\nTime elapsed for Original : | %f \n\n', oend);

NSS = cell(1,9);
Seori = cell(1,2);
Seori(1) = {full(diag(Aori))};
Seori(2) = {full(diag(A0ori))};
NSS(1) = {Seori};
NSS(2) = {oend};
NSS(3) = {i};
NSS(4) = {neg};
NSS(5) = {eps < meps};
NSS(6) = {snfaces};
NSS(7) = {stopfile};
NSS(8) = {Acond};
NSS(9) = {X};

Sori = NSS;

scalemeps = 1;
tmeps = meps * scalemeps^scale;

%% Second, run the algorithm to all scales - MS
msstart = tic;
J=[];
SS = cell(1,finer+1); %preallocate the cell array
for o=mino:finer
    msostart = tic;
    
    tmaxsnfaces = maxsnfaces;
    
    % Lead FIeld
    G = ico(o).cfg.rlf;
    NJs = size(G,2);

    % Initial \alpha
    % If already in the second iteration, use the previous \alpha
    % to initialize this iteration - Multi-Scale Algorithm
    A = ta0*speye(NJs, NJs);
    Aold = A;
    
    if (o > mino)
        tw = tic;

        % Initialize ALpha with values from last iteration and contaminate
        % neighbors with decreased values
%         PrevAI = 1./full(diag(PrevA));
%         Fmap = ico(o).wh.facesmap;
%         for p = 1:NJs
%             Aold(p,p) = PrevAI(Fmap(p));
% 
%             v = ico(o).wh.vizinhos(p,:);
%             nv = size(v,2);
%             for q=1:nv
%                 tnv = size(ico(o).wh.vizinhos(v(q), :),2);
%                 A(v(q),(q)) = A(v(q),v(q)) + Aold(p,p)/tnv;
%             end
%             
%             v2 = cell2mat(ico(o).wh.vizinhos2(p,:));
%             nv2 = size(v2,1);
%             for q=1:nv2
%                 tnv2 = size(cell2mat(ico(o).wh.vizinhos2(v2(q),:)),1);
%                 A(v2(q),(q)) = A(v2(q),v(q)) + Aold(p,p)/tnv2;
%             end
%         end
%         A = 1 ./ A;
%         Aold = 1 ./ Aold;

        % First, contaminate near neigbohrs in the previous dimention
        PrevAI = 1./full(diag(PrevA));
        PrevAIold = PrevAI;
        for p = 1:size(PrevAI,1)
            
            v = ico(o-1).wh.vizinhos(p,:);
            nv = size(v,2);
            for q=1:nv
                tnv = size(ico(o-1).wh.vizinhos(v(q), :),2);
                PrevAI(v(q)) = PrevAI(v(q)) + PrevAIold(p)/tnv;
            end
            
            v2 = cell2mat(ico(o-1).wh.vizinhos2(p,:));
            nv2 = size(v2,1);
            for q=1:nv2
                tnv2 = size(cell2mat(ico(o-1).wh.vizinhos2(v2(q),:)),1);
                PrevAI(v2(q)) = PrevAI(v2(q)) + PrevAIold(v2(q))/tnv2;
            end
        end
        
        % Now, copy this values to the faces spawned in this new dimension
        PrevA = 1./PrevAI;
        PrevAold = 1./PrevAIold;
        Fmap = ico(o).wh.facesmap;
        for p = 1:NJs
            Aold(p,p) = PrevAold(Fmap(p));
            A(p,p) = PrevA(Fmap(p));
        end        
        
        wt = toc(tw);
        fprintf('Tempo para inicializar alphas de %d : %f\n', o, wt);
        
    end

    A0 = A;

    % estimate the noise
%    [~, ~, SGini] = renorm_VB_sparse8(X, G, A, A0, g0, SG);
    
    % Compute J and alpha step until the free energy converges
    if DEBUG
        fAid = fopen([batchdir filesep 'S' num2str(o) '.dat'], 'w');
        fJid = fopen([batchdir filesep 'J' num2str(o) '.dat'], 'w');
    end
    
    prevmalpha = 0;
    prevAcond = 0;
    betams = 0;
    J=[];
    neg=false;
    snfaces = 0;
    prevnfaces = 0;
    stopfile=0;
    for i=1:n
        
        % test stop file
        stopfile=fopen('stopinv');
        if (stopfile > 0)
            warning('Stoping because of stopfile!'); %#ok<WNTAG>
            break;
        end

        % write debug
        if DEBUG
            fprintf(fAid, '%4.16f ', full(diag(A)));
            fprintf(fAid, '\n');
        end
        
        % Save this state
        AP = A;
        JP = J;
        betamsP = betams;
        
        % Try new iteration
        [J, A, betams] = renorm_VB_sparse8(X, G, A, A0, g0, SG);

        % write debug
        if DEBUG
            fprintf(fJid, '%4.16f ', full(J(:,HT)));
            fprintf(fJid, '\n');
        end
        
        DA = full(diag(A));
%         [~, av] = hist(DA, 100);
%         malpha = mean(DA(DA > av(alphalim)));
%        malpha = mean(DA);
        ga = sort(DA);
        malpha = mean(ga(1:ceil(NJs*0.1)));
        maxalpha = max(DA);
        minalpha = min(DA);
        
        Acond = abs(log10(condest(A)));
        Acond1mom = abs(Acond - prevAcond);
        
        fprintf('%d step %d (%1.1e %1.1e)\t', NJs, i, ta0, tg0);
        fprintf('Alpha : > %e %e %e < | ', malpha, maxalpha, minalpha);
        fprintf('Beta : %e %e \n', betams, abs(betams-betamsP));
                
        % Check for numerical errors
        if (minalpha < 0)
            A = AP;
            J = JP;
            neg = 1;
            fprintf('WARNING! alpha or lambda assumed negative values! Keeping previous calculation\n');
            break;
        end
        
        % Check for stop criterea
        aeps = abs(prevmalpha-malpha);
        eps = aeps/abs(malpha);
        [jm, ~] = hist(abs(mean(J,2)), treshold);
        nfaces = sum(jm(histcut:end),2);
        if (nfaces == prevnfaces)
            snfaces = snfaces+1;
        else
            snfaces = 0;
        end

        maxJ = full(mean(max(abs(J))));
        meanJ = full(mean(mean(abs(J))));

%        scalequit = nfaces/NJs;
        fprintf ('Alpha Change: %e (%1.0e) |', eps, tmeps);
        fprintf (' J (%f) : %d (%d) {%d} | [%f (%f)] | ', 1-histcut/treshold, nfaces, snfaces, ceil(tmaxsnfaces), maxJ, meanJ);
        fprintf ('Acond : {%f (15) (%e)} \n', Acond, Acond1mom);
        
        oki = true;
        if o == finer
            oki = i > 100;
        end
        
        if (prevmalpha > 0 && oki && (eps < tmeps) )% || ...
%                (o < finer && minalpha < minalphaori) ...
            % || Acond > 15)% || (scalequit < 0.1 && o < finer))
            break;
        end
        
        prevmalpha = malpha;
        prevnfaces = nfaces;
        prevAcond = Acond;

    end
    if (i == n)
        fprintf('Multi-Scale iteration %d exit for maximum iteration. Considering change minimum alpha variation\n', o);
    end
    fprintf ('Fim da escala %d\n\n', o);
    if DEBUG
        fclose(fAid);
        fclose(fJid);
    end
    
    msoend = toc(msostart);
    NS = cell(1,8);
    SeSini = cell(1,3);
    SeSini(1) = {full(diag(A))};
    SeSini(2) = {full(diag(Aold))};
    SeSini(3) = {full(diag(A0))};
    NS(1) = {SeSini};
    NS(2) = {msoend};
    NS(3) = {i};
    NS(4) = {neg};
    NS(5) = {eps < tmeps};
    NS(6) = {snfaces};
    NS(7) = {Acond};
    NS(8) = {stopfile};
    
    SS(o) = {NS};
    
    PrevA = A;
    
    tmeps = tmeps/scalemeps;
    ta0 = ta0/scalea;
    tg0 = tg0/scaleg;
    
end

msend = toc(msstart);
fprintf('\nTime elapsed for Multi-Scale: | %f \n\n', msend);

SS(finer+1) = {msend};


% compare with the default MNE
Amne = a0*speye(NJs, NJs);
G = sparse(G);
M1 = Amne * G';
M2 = ( G * Amne * G' + SG);
M3 = M1 / M2;
Jmmq = M3 * X;
end
