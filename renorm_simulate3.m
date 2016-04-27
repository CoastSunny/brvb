function [cfg radius]= renorm_simulate3(vertices, faces, ico, varargin)


go = true;
if nargin > 3
    go = varargin{1};
end

limsup=1000;
if nargin > 4
    limsup = varargin{2};
end

liminf=5;
if nargin > 5
    liminf = varargin{3};
end

smear = true;
if nargin > 6
    smear = varargin{4};
end

cfg = ico.cfg;
ndips = size(cfg.dip.mom,1);

% extended version
emom = cfg.dip.mom;
cfg.dip.j = zeros(ndips,1);

% find two random dipoles inside de minimum range (or 1000 if not
% specified)
if limsup > liminf
    % simulate 2 random dipoles
    active = 2;
    [~, hot radius] = renorm_rand_dip(cfg, limsup, liminf, vertices, ico, go);
else
    % simulate 1 random dipole
    active = 1;
    one = ceil(ndips * rand(1,1)); % one random dipole
    hot = one;
end

cfg.dip.j(hot) = sign(rand(active,1)-0.5) .* (1 - 0.7 * rand(active,1));
cfg.dip.hot = hot;

if smear
    for h=1:size(hot,2)
        v1 = ico.wh.vizinhos(hot(h),:);
        tnv1 = size(v1,2);
        v2 = cell2mat(ico.wh.vizinhos2(hot(h),:));
        tnv2 = size(v2,1);
        for q=1:tnv1
            cfg.dip.j(v1(q)) = cfg.dip.j(v1(q)) + cfg.dip.j(hot(h))/tnv1;
        end

        for q=1:tnv2
            cfg.dip.j(v2(q)) = cfg.dip.j(v2(q)) + cfg.dip.j(hot(h))/tnv2;
        end

    end

end

% PROPAGATE IN TIME
T = 51;
nj = zeros(ndips, T);
x =linspace(0,pi(),T);
for i=1:T
     nj(:,i) = sin(x(i)) * cfg.dip.j;
end
cfg.dip.j = nj;

% extended version
x = size(emom,1);
y = size(emom,2);
emom = reshape(emom', x * y, 1);
emom = repmat(emom,1,T);
for i=1:ndips
    emom(i*y,:) = emom(i*y,:) .* cfg.dip.j(i,:);
    emom(i*y-1,:) = emom(i*y-1,:) .* cfg.dip.j(i,:);
    emom(i*y-2,:) = emom(i*y-2,:) .* cfg.dip.j(i,:);
%   emom(i,:) = emom(i,:) .* cfg.dip.j(i);
end
cfg.dip.emom = emom;

% reduced version
cfg.dip.rmom = cfg.dip.j;

%compute the eletric potential for the selected dipoles
cfg.pot = cfg.lf * cfg.dip.emom;
cfg.rpot = cfg.rlf * cfg.dip.rmom;
cfg.diff = cfg.pot - cfg.rpot;

% plot the 3-D distribution of the potential over the sphere surface ?
if go
    init = ceil (T/2);
    
    cfg.elec.tri = convhulln(cfg.elec.pnt);
    f2 = figure;

    patch(struct('vertices', vertices, 'faces', faces, 'FaceColor', 'g', 'FaceAlpha', 0, 'EdgeAlpha', 0.1));
    hold all;
    patch('faces', cfg.elec.tri, 'vertices', cfg.elec.pnt, 'FaceVertexCData', cfg.rpot(:,init), 'FaceColor', 'interp', 'FaceAlpha', 0.5);

%    patch(struct('vertices', vertices, 'faces', faces(hot), 'FaceColor', 'g', 'FaceAlpha', 1));
    
    for i=1:size(faces,1)
         if ismember(i, hot)
            fprintf('imprime seta...\n');
            v = cfg.dip.pos(i,:) + 50*cfg.dip.rmom(i)*cfg.dip.mom(i,:);
            vectarrow(cfg.dip.pos(i,:), v(1,:));
        end
    end

    h = gca;
    axis (h, 'equal');
    axis (h,'vis3d');
    cameramenu (f2);

% 
%     f1 = figure;
% 
%     patch(struct('vertices', vertices, 'faces', faces, 'FaceColor', 'g', 'FaceAlpha', 0.5));
%     hold all;
%     patch('faces', cfg.elec.tri, 'vertices', cfg.elec.pnt, 'FaceVertexCData', cfg.pot(:,init), 'FaceColor', 'interp', 'FaceAlpha', 0.5);
% 
%     for i=1:size(faces,1)
%         if ismember(i, hot)
%             v = cfg.dip.pos(i,:) + 50*cfg.dip.rmom(i)*cfg.dip.mom(i,:);
%             vectarrow(cfg.dip.pos(i,:), v(1,:));
%         end
%     end
% 
%     h = gca;
%     axis (h, 'equal'); 
%     axis (h, 'vis3d'); 
%     cameramenu (f1);

    f3 = figure;

%     for j=1:T
        clf(f3);

        patch(struct('vertices', vertices, 'faces', faces, 'FaceColor', 'g'));
        hold all;
        patch('faces', cfg.elec.tri, 'vertices', cfg.elec.pnt, 'FaceVertexCData', cfg.diff(:,ceil(T/2)), 'FaceColor', 'interp', 'FaceAlpha', 0.5);
 
        drawnow;

        t=0.1;
%         if (nargin > 3)
%             t = varargin{1};
%         end

        pause(t);

%     end

    h = gca;
    axis (h, 'equal'); 
    axis (h, 'vis3d'); 
    cameramenu(f3);
    
end
