function [STATE, j] = renorm_load3(datapath, file, alpha)

STATE = [];
j=0;
if str2num(file.name)
    j = str2num(file.name);
    path = [datapath file.name '/STATE.mat'];
    S = {open(path)};
    
    if isfield(S{1}.STATE.cfg.dip, 'hot')
        STATE.cfg.dip.hot = S{1}.STATE.cfg.dip.hot;
    else
        STATE.cfg.dip.hot = find_hot(S{1}.STATE.cfg.dip.rmom);
    end
    STATE.cfg.rpot = S{1}.STATE.cfg.rpot;
    %        size(S{1}.STATE.cfg.dip.rmom)
    STATE.cfg.dip.rmom = S{1}.STATE.cfg.dip.rmom;
%     STATE.cfg.dip.emom = S{1}.STATE.cfg.dip.emom;
%     STATE.cfg.dip.j = S{1}.STATE.cfg.dip.j;
    
    STATE.dip = S{1}.STATE.dip;
    STATE.dipori = S{1}.STATE.dipori;
    STATE.dipmmq = S{1}.STATE.dipmmq;
    if isfield(S{1}.STATE, 'X')
        STATE.X = S{1}.STATE.X;
    end
    
    if (alpha)
        for k=1:5
            if (~isempty(S{1}.STATE.SS{k}))
                STATE.SS{k}{1}{3} = S{1}.STATE.SS{k}{1}{3};
                STATE.SS{k}{1}{1} = S{1}.STATE.SS{k}{1}{1};
            end
        end
        STATE.SSOri{1}{1} = S{1}.STATE.SSOri{1}{1};
    end
    %
    %
    %         STATE.cfg.dip.pos = S.STATE.cfg.dip.pos;
    %         STATE.SSOri{2} = S.STATE.SSOri{2};
    %         order = size(S.STATE.SS,1);
    %         STATE.SS{order+1} = S.STATE.SS{order+1};
end

end

function h = find_hot(z)
    t = round(size(z,2)/2);
    x = z(find(z(:,t)),t);
    y = unique(x);
    a = histc(x, y);
    hv = y(a==1);
    if length(hv) > 2, error ('Can''t find original dipoles'); end
    h = zeros(2,1);
    
%     figure 
%     hold on
%     c = find(z(:,t));
%     for ic=1:length(c)
%         plot(z(c(ic),:), 'r')
%     end

    for id = 1:2
        h(id) = find(z(:,t)==hv(id)); 
%         plot(z(h(id),:))
%         hold on
    end
    
end