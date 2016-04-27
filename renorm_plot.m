function figures = renorm_plot(surface, ico, STATE, varargin)

fprintf('Teste\n');

histcut = 2;
if nargin > 3
    histcut = varargin{1};
%    order = varargin{1};
end

inv = true;
if nargin > 4
    inv = varargin{2};
%    order = varargin{1};
end

order = 5;
if nargin > 5
    order = varargin{3};
end

palpha = true;
if nargin > 6
    palpha = varargin{4};
end

if isfield(STATE, 'STATE')
    STATE = STATE.STATE;
end
  
real = STATE.cfg.dip.rmom;
dipori = STATE.dipori;
dipmmq = STATE.dipmmq;
dip = STATE.dip;

global EdgeAlpha;
EdgeAlpha = 0.1;
global FaceAlpha;
FaceAlpha = 1;

global totalndips;
totalndips = size(real,1)
S.totalndips = totalndips;
global totaltime;
totaltime = size(real,2)
S.totaltime = totaltime;

global init;
init = 1;
if totaltime > 1
    init = ceil(totaltime/2);
end



realnz = find(real(:, init));
ndips = size(realnz,1);

[~, msdip, ~, ~, ~, treshold] = renorm_set_treshold(ndips, dip(:,init), histcut); 
[~, vbdip] = renorm_set_treshold(ndips, dipori(:,init), histcut);
[~, mmqdip] = renorm_set_treshold(ndips, dipmmq(:,init), histcut);

%%%% Real
alpha = FaceAlpha * zeros(totalndips,1);
alpha(realnz) = 1;

f=0;
f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('Real', surface, ico, order, real, alpha, histcut, treshold, size(realnz,1), ndips);

%%%% MS
alpha = FaceAlpha * zeros(totalndips,1);
alpha(msdip) = 1;
size(alpha)
f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('MG', surface, ico, order, dip, alpha, histcut, treshold, size(msdip,1), ndips);

%%%% VB
alpha = FaceAlpha * zeros(totalndips,1);
alpha(vbdip) = 1;

f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('VB', surface, ico, order, dipori, alpha, histcut, treshold, size(vbdip,1), ndips);

%%%% MNE

alpha = FaceAlpha * zeros(totalndips,1);
alpha(mmqdip) = 1;

f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('MNE', surface, ico, order, dipmmq, alpha, histcut, treshold, size(mmqdip,1), ndips);

%%%%% coarse grained 

if palpha
    for i=1:order
        
        if (isfield(STATE, 'SS') && ~isempty(STATE.SS{i}))
            
            nhistcut = histcut;

            X = STATE.SS{i}{1}{3};
            if inv
                X = 1./X;
            end
            
            n = size(X,1);
            totalndips = n;
            totaltime = 1;
            init = 1;
            
            alpha = FaceAlpha * zeros(n,1);
            [~, vba, ~, ~, ~, treshold] = renorm_set_treshold(ndips, X, nhistcut);
%             if i == 1
%                 nhistcut = 0;
%                 alpha = 1;
%                 X = 1;
%             else
                alpha(vba) = 1;
%             end
            
            f = f + 1; [figures(f,1), figures(f,2)] = cria_patch(['Ordem MS ' num2str(i) ' Inicial'], surface, ico, i, X, alpha, nhistcut, treshold, size(vba,1), ndips);
            
            nhistcut = histcut;
%             if i == order
%                 nhistcut = 181;
%             end
            
            X = STATE.SS{i}{1}{1};
            if inv
                X = 1./X;
            end
            [~, vba, ~] = renorm_set_treshold(ndips, X, nhistcut);
            alpha = FaceAlpha * zeros(n,1);
            alpha(vba) = 1;
            
            f = f + 1; [figures(f,1), figures(f,2)] = cria_patch(['Ordem MS ' num2str(i) ' Final'], surface, ico, i, X, alpha, nhistcut, treshold, size(vba,1), ndips);
            
        end
    end
    
    if isfield(STATE, 'SSOri')
        if (iscell(STATE.SSOri{1}))
            X = STATE.SSOri{1}{1};
        else
            X = diag(STATE.SSOri{1});
        end
        
        if inv
            X = 1./X;
        end
        
        n = size(X,1);
        totalndips = n;
        totaltime = 1;
        init = 1;
        
        f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('1/Alpha VB', surface, ico, i, ones(1,size(X,2)), ones(1,size(X,2)), histcut, treshold, size(vba,1), ndips);

%         histcut = 60;
        
        [~, vba, ~] = renorm_set_treshold(ndips, X, histcut);
        alpha = FaceAlpha * zeros(n,1);
        alpha(vba) = 1;
        
        f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('1/Alpha VB', surface, ico, i, X, alpha, histcut, treshold, size(vba,1), ndips);
        
    end
    
end

end

function [fig, p] = cria_patch(name, surface, ico, order, X, alpha, histcut, treshold, eda, ndips)

    S.dip = X;

    S.ndips = ndips;
    slidebegin = 20;
    slidesize = 500;

    global EdgeAlpha;
    global FaceAlpha;
    S.FaceAlpha = FaceAlpha;

    global totalndips;
    S.totalndips = totalndips;
    global totaltime;
    S.totaltime = totaltime;

    global init;
    
    fig = figure;
    S.f = fig;
    p = patch('vertices', surface.vertices, 'faces', ico(order).wh.faces, 'FaceVertexCData', X(:,init), ...
        'FaceColor', 'flat', 'DisplayName', name, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha,...
        'EdgeAlpha', EdgeAlpha);
    S.p = p;
    legend show;
    axis off;
    axis equal;
    colorbar;
    if min(X(:,init)) < max(X(:,init))
        caxis([min(X(:,init)) max(X(:,init))]);
    end
    cameramenu;
    % transparency
    S.sla = uicontrol('style','slide',...
                     'unit','pix',...
                     'position',[slidebegin 10 slidesize 30],...
                     'min',1,'max',1000,'val', histcut);
    S.edh = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[slidesize+slidebegin+10 10 150 30],...
                 'fontsize',16,...
                 'string',num2str(0.1*(treshold - histcut)));
    S.eda = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[slidesize+slidebegin+10+150+10 10 150 30],...
                 'fontsize',16,...
                 'string',num2str(eda));

    if totaltime > 1
        % actual current
        S.slt = uicontrol('style','slide',...
                        'unit','pix',...
                        'position',[slidebegin 40 slidesize 30],...
                        'min',1,'max',totaltime,'val', init);

        S.edt = uicontrol('style','edit',...
                     'unit','pix',...
                     'position',[slidesize+slidebegin+10 40 150 30],...
                     'fontsize',16,...
                     'string', num2str(init));
        S.edmc = uicontrol('style','edit',...
                     'unit','pix',...
                     'position',[slidesize+slidebegin+10+150+10 40 150 30],...
                     'fontsize',16,...
                     'string', num2str(max(X(:,init))));
        set(S.slt,'call',{@time_call,S});
    else
        S.slt = 0;
        S.edt = 0;
        S.edmc = 0;
    end
             
    set(S.sla,'call',{@alpha_call,S});

%     fullscreen = get(0,'ScreenSize');
%     set(S.f, 'Position',[0 -50 fullscreen(3) fullscreen(4)]);
end


function [] = alpha_call(varargin)
% Callback for the edit box and slider.

[h,S] = varargin{[1,3]};  % Get calling handle and structure.

histcut = ceil(get(h,'value'));
if S.slt
    t = ceil(get(S.slt,'value'));
else
    t = 1;
end

[~, dip, ~, ~, ~, treshold] = renorm_set_treshold(S.ndips, S.dip(:,t), histcut);
alpha = S.FaceAlpha * zeros(S.totalndips,1);
alpha(dip) = 1;

if (isempty(find(alpha==0)))
    set(S.p, 'FaceAlpha', 1);
else
    set(S.p, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha);
end


set(S.edh,'string',num2str(0.1*(treshold - histcut)));
set(S.eda,'string',num2str(size(dip,1)));
end
 
 
function [] = time_call(varargin)
% Callback for the edit box and slider.

[h,S] = varargin{[1,3]};  % Get calling handle and structure.

t = ceil(get(h,'value'));
set(S.p, 'FaceVertexCData', S.dip(:,t));

if (S.edmc)
    set(S.edmc,'string',num2str(max(S.dip(:,t))));
end
if (S.edt)
    set(S.edt,'string',num2str(t));
end

alpha_call(S.sla, 0, S);

end
