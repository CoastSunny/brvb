function [p, ss, jj] = renorm_plot_alpha(wh, ico, o, ss, jj, varargin)

ori=0;
if (nargin > 7)
    ori = varargin{3};
end
if (ischar(ss))
    d = ss;
    if ori
        ss = dlmread([d filesep 'Sori.dat'], ' ');
        jj = dlmread([d filesep 'Jori.dat'], ' ');
    else
        ss = dlmread([d filesep 'S' num2str(o) '.dat'], ' ');
        jj = dlmread([d filesep 'J' num2str(o) '.dat'], ' ');
    end
end

histcut = 20;
EdgeAlpha = 0.1;
FaceAlpha = 0.5;

fig = figure;
maxfigsize(fig);
ns = size(ss,1);

max = 1;
wfaces = [];
if (nargin > 5)
    max = varargin{1};
end

inv = true;
if (nargin > 6)
    inv = varargin{2};
end

k=1;
for i = 1:ns
    
    if (rem(i,max) == 0 || i == 1)
        
        l = 1;
        if (max+i>ns)
            l=0;
        end
        
        k=k+1;
        clf(fig);
        
        X = ss(i,:)';
        X(end) = [];
        
%         n = size(X,1);
%         [~, vba, ~] = renorm_set_treshold(n, X, histcut);
%         alpha = FaceAlpha * ones(n,1);
%         alpha(vba) = 1;
        
        %     X=ones(n,1);
        %     X(find(STATE.SS{i}{1}{2} > 1.0e-8)) = 1; %#ok<*FNDSB>
        %     X(find(STATE.SS{i}{1}{2} < 1.0e-8)) = 0;
        
        %        'FaceColor', 'flat', 'DisplayName', ['Ordem ' num2str(o) ' Inicial'], 'FaceAlpha', 0.5, ...

%         [~, vba, ~] = renorm_set_treshold(ndips, X, histcut);
%         alpha = FaceAlpha * ones(n,1);
%         alpha(vba) = 1;
        if isnumeric(o)
            faces = ico(o).wh.faces;
        else
            faces = ico(5).wh.faces;
        end

        if inv
            X = 1./X;
        end
       
        h1 = subplot(1,2,1);
        p = patch('vertices', wh.orig.vertices, 'faces', faces, 'FaceVertexCData', X, ...
            'FaceColor', 'flat', 'DisplayName', ['S' num2str(o) ' batch'], ... %, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha, ...
            'EdgeAlpha', EdgeAlpha);
        %    patch('vertices', wh.orig.vertices, 'faces', ico(o).wh.faces(wfaces,:), 'FaceColor', 'k', 'FaceAlpha', 1);
        if l && ~(nargin > 8 && i ==1)
            rotate(p, [1 0 0], -90);
            rotate(p, [0 1 0], -10*k);
        end
        
        legend show;
        axis off;
        axis equal;
        colorbar;
        cameramenu;
%        caxis([min(X), max(X)]);
        
        Y = jj(i,:)';
        Y(end) = [];

        h2 = subplot(1,2,2);
        p = patch('vertices', wh.orig.vertices, 'faces', faces, 'FaceVertexCData', Y, ...
            'FaceColor', 'flat', 'DisplayName', ['J' num2str(o) ' batch'], ... %, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha, ...
            'EdgeAlpha', EdgeAlpha);
        %    patch('vertices', wh.orig.vertices, 'faces', ico(o).wh.faces(wfaces,:), 'FaceColor', 'k', 'FaceAlpha', 1);
        if l && ~(nargin > 8 && i ==1)
            rotate(p, [1 0 0], -90);
            rotate(p, [0 1 0], -10*k);
        end
        
        legend show;
        axis off;
        axis equal;
        colorbar;
        cameramenu;
        
        drawnow;
        if (nargin > 8 && i ==1)
            drawnow;
            pause(varargin{4})
%             %    else
%             %        pause(0.01);
        end
        
    fprintf('%d de %d\n', i, ns);
    end
    
end

cameramenu;
