function [index hot radius] = renorm_rand_dip(cfg, limsup, liminf, vertices, ico, go)

if nargin < 6, go = false; end

totalndips = size(cfg.dip.pos,1);
index = [];
while (size(find(index),1) < 2)
    index = false(totalndips,1);
    one = ceil(totalndips * rand(1,1)); % two random dipoles
    index(one) = 1;
    radius = ((limsup-liminf) * rand(1,1))+liminf;
    snext = randblock(linspace (1, totalndips, totalndips),1);
    for j=1:totalndips
        two = snext(j);
        index(two) = 1;
        d = sqrt(sum(diff(cfg.dip.pos(index,:)).^2));
        if (d+1 > radius && d-1 < radius && snext(j) ~= one)
            fprintf ('vai sair. one: %d two : %d\n', one, two);
            break;
        else
            index(two) = 0;
        end
    end
end

fprintf ('Simulation results : radius %2.1f d %2.1f\n', radius, d)
hot = [one two];

if go
    
    EdgeAlpha = 0.1;
    FaceAlpha = 0.5;
    
    alpha = FaceAlpha * ones(totalndips,1);
    alpha(index) = 1;
    x = zeros(size(cfg.dip.pos,1),1);
    x(index)=1;
    x(not(index))=0;
    
    figure;
    patch('vertices', vertices, 'faces', ico.wh.faces, 'FaceVertexCData', x, ...
        'FaceColor', 'flat', 'DisplayName', 'Random Dips', 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alpha,...
        'EdgeAlpha', EdgeAlpha)
    legend show;
    axis off;
    axis equal;
    colorbar;
    caxis([-1. 1]);
    cameramenu;
end