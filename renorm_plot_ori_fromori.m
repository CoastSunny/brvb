 function figures = renorm_plot_ori_fromori(surface, ico, STATE, varargin)


global surf;
surf= 1; % 1 = original, 2 = inflated, 3 = sphere

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

normalize = true;
if nargin > 7
    normalize = varargin{5};
end

if isfield(STATE, 'STATE')
    STATE = STATE.STATE;
end
  
realori = STATE.cfg.dip.rmom;
real = zeros(size(ico(order).wh.faces,1), size(STATE.cfg.dip.rmom,2));
for j=1:size(STATE.cfg.dip.rmom,2)
    realorifaces = find(realori(:,j));
    if ~isempty(realorifaces)
        real(ico(order).wh.facesmapori(realorifaces(1)),j) = realori(realorifaces(1),j);
        real(ico(order).wh.facesmapori(realorifaces(2)),j) = realori(realorifaces(2),j);
    end
end

dipori = STATE.dipori;
dipmmq = STATE.dipmmq;
dip = STATE.dip;

global EdgeAlpha;
EdgeAlpha = 0;
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
alpha = zeros(totalndips,1);
alpha(realnz) = FaceAlpha;

f=0;
f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('Real', surface, ico, order, real, alpha, histcut, treshold, size(realnz,1), ndips, normalize);

%%%% MS
alpha = zeros(totalndips,1);
alpha(msdip) = FaceAlpha;
size(alpha)
f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('MG', surface, ico, order, dip, alpha, histcut, treshold, size(msdip,1), ndips, normalize);

%%%% VB
alpha = zeros(totalndips,1);
alpha(vbdip) = FaceAlpha;

f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('VB', surface, ico, order, dipori, alpha, histcut, treshold, size(vbdip,1), ndips, normalize);

%%%% MNE

alpha = zeros(totalndips,1);
alpha(mmqdip) = FaceAlpha;

f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('MNE', surface, ico, order, dipmmq, alpha, histcut, treshold, size(mmqdip,1), ndips, normalize);

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
            
            alpha = zeros(n,1);
            [~, vba, ~, ~, ~, treshold] = renorm_set_treshold(ndips, X, nhistcut);
%             if i == 1
%                 nhistcut = 0;
%                 alpha = 1;
%                 X = 1;
%             else
                alpha(vba) = FaceAlpha;
%             end
            
            f = f + 1; [figures(f,1), figures(f,2)] = cria_patch(['Ordem MS ' num2str(i) ' Inicial'], surface, ico, i, X, alpha, nhistcut, treshold, size(vba,1), ndips, normalize);
            
            nhistcut = histcut;
%             if i == order
%                 nhistcut = 181;
%             end
            
            X = STATE.SS{i}{1}{1};
            if inv
                X = 1./X;
            end
            [~, vba, ~] = renorm_set_treshold(ndips, X, nhistcut);
            alpha = zeros(n,1);
            alpha(vba) = FaceAlpha;
            
            f = f + 1; [figures(f,1), figures(f,2)] = cria_patch(['Ordem MS ' num2str(i) ' Final'], surface, ico, i, X, alpha, nhistcut, treshold, size(vba,1), ndips, normalize);
            
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
        
        [~, vba, ~] = renorm_set_treshold(ndips, X, histcut);
        
        f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('1/Alpha VB', surface, ico, order, FaceAlpha * ones(size(X,1), 1), ones(size(X,1),1), histcut, treshold, size(vba,1), ndips, normalize);

%         histcut = 60;
        

        alpha = zeros(n,1);
        alpha(vba) = FaceAlpha;
        
        f = f + 1; [figures(f,1), figures(f,2)] = cria_patch('1/Alpha VB', surface, ico, order, X, alpha, histcut, treshold, size(vba,1), ndips, normalize);
        
    end
    
end

end

function [fig, p2] = cria_patch(name, surface, ico, order, X, alpha, histcut, treshold, eda, ndips, normalize)

    S.dip = X;

    S.ndips = ndips;
    slidebegin = 20;
    slidesize = 500;

    global surf;
    global EdgeAlpha;
    global FaceAlpha;
    S.FaceAlpha = FaceAlpha;

    global totalndips;
    S.totalndips = totalndips;
    global totaltime;
    S.totaltime = totaltime;

    global init;
    
    % normalize X
    if normalize
        nX = max(max(abs(X)));
        X = X ./ nX;
    end
    
    if size(X,2) > 1
        tX = X(:,init);
    else
        tX = X(:,1);
    end
    
    fig = figure;
    set(fig, 'Visible', 'off')
    
    S.f = fig;
    S.smap = ico(order).smap;
    xcolor = zeros(size(ico(order).wh.facesmapori));
    xalpha = zeros(size(ico(order).wh.facesmapori));
    
    for i=1:totalndips
        xcolor(S.smap{i}) = tX(i);
%         if mod(i,5)
            xalpha(S.smap{i}) = alpha(i);
%         end
    end
    S.xcolor = xcolor;
    S.xalpha = xalpha;
    
    switch surf
        case 1
            s = surface.orig;
        case 2
            s = surface.inflated;
        case 3
            s = surface.sphere;
    end
    c = surface.curv.cdata;

    % CURV COLOR CODE 1
    c = c - min(c);
    c = 0.5 * c ./ max(c);

    % CURV COLOR CODE 2
%     c = c - min(c);
%     c = c ./ max(c);
%     mc = mean(c);
%     c(c>mc) = 0.6;
%     c(c<mc) = 0.4;
     
    hold off;
    if surf == 2
        p1 = patch('vertices', s.vertices, 'faces', s.faces, 'FaceVertexCData', [c c c], ...
            'FaceColor', 'flat', 'DisplayName', name, 'EdgeAlpha', EdgeAlpha); %..., 'FaceAlpha', 'flat', 'FaceVertexAlphaData', 0.1);
        hold on;
        p2 = patch('vertices', s.vertices, 'faces', s.faces, 'FaceVertexCData', xcolor, ...
            'FaceColor', 'flat', 'DisplayName', name, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', xalpha,...
            'EdgeAlpha', EdgeAlpha);
        S.p1 = p1;
        S.p2 = p2;
    else
        p2 = patch('vertices', s.vertices, 'faces', s.faces, 'FaceVertexCData', xcolor, ...
            'FaceColor', 'flat', 'DisplayName', name, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', xalpha,...
            'EdgeAlpha', 0.1);
        S.p2 = p2;
    end
    saveas(gcf, [name '.fig'])
    S.surf = surf;
%     legend show;
    axis off;
    axis equal;
    colorbar;
    if min(tX) < max(tX)
        caxis([min(tX) max(tX)]);
    end
    cameramenu;
    
%     rotate(p1, [0 1 0], -90)
    rotate(p2, [0 1 0], -90)

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
                     'string', num2str(max(tX)));
        set(S.slt,'call',{@time_call,S});
    else
        S.slt = 0;
        S.edt = 0;
        S.edmc = 0;
    end
             
    set(S.sla,'call',{@alpha_call,S});

    fullscreen = get(0,'ScreenSize');
    set(S.f, 'Position',[0 -50 fullscreen(3) fullscreen(4)]);
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

alpha = zeros(size(S.xalpha));
alpha(cell2mat(S.smap(dip))) = S.FaceAlpha;


if (isempty(find(alpha==0)))
    set(S.p2, 'FaceAlpha', S.FaceAlpha);
else
    set(S.p2, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha);
end


set(S.edh,'string',num2str(0.1*(treshold - histcut)));
set(S.eda,'string',num2str(size(dip,1)));
end
 
 
function [] = time_call(varargin)
% Callback for the edit box and slider.

[h,S] = varargin{[1,3]};  % Get calling handle and structure.

t = ceil(get(h,'value'));

xcolor = zeros(size(S.xcolor));
for i=1:length(S.dip)
    xcolor(S.smap{i}) = S.dip(i,t);
end
    
set(S.p2, 'FaceVertexCData', xcolor);

if (S.edmc)
    set(S.edmc,'string',num2str(max(S.dip(:,t))));
end
if (S.edt)
    set(S.edt,'string',num2str(t));
end

% alpha_call(S.sla, 0, S);

end
