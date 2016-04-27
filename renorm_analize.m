function S = renorm_analize(MVBdata, ico, states, o)

R = 3;

if (ischar(states))
    files = dir(states);
    n = size(files,1)-2;
else
    n = size(states,2);
end

ms = zeros(n,R+1);
vb = zeros(n,R+1);
mne = zeros(n,R+1);
jr = zeros(n,R+1);

msn = zeros(n,R+1);
vbn = zeros(n,R+1);
mnen = zeros(n,R+1);
jrn = zeros(n,R+1);

x = zeros(1, R+1);

dms = zeros(n,3);
dvb = zeros(n,3);
dmne = zeros(n,3);

errms = zeros(n,R+1);
errvb = zeros(n,R+1);
errmne = zeros(n,R+1);

k=1;

T = 0;
totalndips = size(ico(o).wh.mean,1);
ondips = totalndips;

% 
% xcolor = ones(10240,3);
% figure; p = patch('vertices', MVBdata.wh.inflated.vertices, 'faces', MVBdata.ico(5).wh.faces, 'FaceVertexCData', xcolor, 'FaceColor', 'flat');
% 

if exist('vizinhos.mat', 'file')
    load('vizinhos.mat');
else
    vizinhos = cell(totalndips,R+1);
    for j=1:totalndips
        fprintf('J = %d\n',j);
        vizinhos(j,1) = {j};
        for r=2:R+1
            v = vizinhos{j,r-1};
            nv = [];
            for rr = 1:size(v,2)
                v1 = ico(o).wh.vizinhos(v(rr),:);
                v1 = reshape(v1,1,size(v1,1)*size(v1,2));
                v2 = cell2mat(ico(o).wh.vizinhos2(v(rr),:))';
                nv = unique([nv v1 v2]);
            end
            vizinhos(j,r) = {nv};
        end
%         xcolor = ones(10240,3);
%         xcolor(cell2mat(vizinhos(j,4)), :) = repmat([0 0 1],length(cell2mat(vizinhos(j,4))), 1);
%         xcolor(cell2mat(vizinhos(j,3)), :) = repmat([0 1 0],length(cell2mat(vizinhos(j,3))), 1);
%         xcolor(cell2mat(vizinhos(j,2)), :) = repmat([1 0 0],length(cell2mat(vizinhos(j,2))), 1);
%         xcolor(cell2mat(vizinhos(j,1)), :) = repmat([1 1 0],length(cell2mat(vizinhos(j,1))), 1);
%         set(p,'FaceVertexCData', xcolor);
    end
    save('vizinhos.mat', 'vizinhos');
end

for i=1:n

    if (ischar(states))
        fprintf('opening %d = %s...\n', i, files(i).name);
        state = renorm_load3(states, files(i), false);
    else
        if (isfield(states{i}, 'STATE'))
            state = states{i}.STATE;
        else
            state = states{i};
        end
    end
    
    if ~isempty(state)
        fprintf('Calculando rho para %d\n', i);

        if T == 0
            dip = mean(state.cfg.dip.rmom,2);
            dipms = mean(state.dip,2);
            dipvb = mean(state.dipori,2);
            dipmne = mean(state.dipmmq,2);
        else
            dip = state.cfg.dip.rmom(:,T);
            dipms = state.dip(:,T);
            dipvb = state.dipori(:,T);
            dipmne = state.dipmmq(:,T);
        end
        
        if size(dip,1) > ondips
            fprintf('Simulation used real surface to generate field. Comparing distance to original dipole instead of it''s placeholder in the %d th surface.\n', o);
            realmeans = MVBdata.cfg.dip.pos;
        else
            realmeans = ico(o).wh.mean;
        end
        totalndips = size(dip,1);
        
        Z = sqrt(sum(dip.^2, 1));
        
        dipn = abs(dip);
        dipn = (dipn-min(dipn))./(max(dipn)-min(dipn));
        Zn = sqrt(sum(dipn.^2, 1));

        dipmsn = abs(dipms);
        dipmsn = (dipmsn-min(dipmsn))./(max(dipmsn)-min(dipmsn));

        dipvbn = abs(dipvb);
        dipvbn = (dipvbn-min(dipvbn))./(max(dipvbn)-min(dipvbn));

        dipmnen = abs(dipmne);
        dipmnen = (dipmnen-min(dipmnen))./(max(dipmnen)-min(dipmnen));

%         errms(k,1) = {dip - dipms};
%         errvb(k,1) = {dip - dipvb};
%         errmne(k,1) = {dip - dipmne};

        real = find(dip);
        ndips = size(real,1);
        [~, realsp] = sort(abs(dip), 'descend');

        % first, calculate the distace
        [~, dipvbsp] = sort(abs(dipvb), 'descend');
        [~, dipmssp] = sort(abs(dipms), 'descend');
        [~, dipmnesp] = sort(abs(dipmne), 'descend');
%         hot = find(ismember(dip, state.cfg.dip.hot));
%         [~, realhotsp] = sort(abs(dip(hot)), 'descend');
        for j=1:2
             dvb(k,j) = sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipvbsp(j),:)).^2));
             dms(k,j) = sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipmssp(j),:)).^2));
             dmne(k,j) = sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipmnesp(j),:)).^2));
        end
        for j=1:ndips
            dvb(k,3) = dvb(k,3) + sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipvbsp(j),:)).^2));
            dms(k,3) = dms(k,3) + sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipmssp(j),:)).^2));
            dmne(k,3) = dmne(k,3) + sqrt(sum((realmeans(realsp(j),:) - ico(o).wh.mean(dipmnesp(j),:)).^2));
        end
        dvb(k,3) = dvb(k,3)/ndips;
        dms(k,3) = dms(k,3)/ndips;
        dmne(k,3) = dmne(k,3)/ndips;

        % now calculate Rho
        xt = zeros(1, R+1);
        
        Zms = zeros(1, R+1);
        Zvb = zeros(1, R+1);
        Zmne = zeros(1, R+1);
        Zr = zeros(1, R+1);
        
        Zmsn = zeros(1, R+1);
        Zvbn = zeros(1, R+1);
        Zmnen = zeros(1, R+1);
        Zrn = zeros(1, R+1);
        
    %     tic;
%         Znb = 0;
%         nnn=0;
%        for j=1:ndips
        if totalndips > ondips
            fprintf('Can''t calculate rho for field generated in the original surface\n');
        else
            for j=1:totalndips

    %            v = realsp(j);
                for r = 1:R+1

                    v = cell2mat(vizinhos(j,r));
                    xt(r) = xt(r) + size(v, 2);

                    % (J_real^i * <J^i>_r) / (Z_real*Z)
                    %J = abs(dip(realsp(j)));
                    J = abs(dip(j));

                    msJr = mean(abs(dipms(v)));
                    vbJr = mean(abs(dipvb(v)));
                    mneJr = mean(abs(dipmne(v)));
                    Jr = mean(abs(dip(v)));

                    Zms(r) = Zms(r) + msJr^2;
                    Zvb(r) = Zvb(r) + vbJr^2;
                    Zmne(r) = Zmne(r) + mneJr^2;
                    Zr(r) = Zr(r) + Jr^2;

                    ms(k,r) = ms(k,r) + msJr * J;
                    vb(k,r) = vb(k,r) + vbJr * J;
                    mne(k,r) = mne(k,r) + mneJr * J;
                    jr(k,r) = jr(k,r) + Jr * J;

                    % (J_real^i * <Jn^i>_r) / (Z_real*Zn)
                    J = abs(dipn(j));

                    msJrn = mean(abs(dipmsn(v)));
                    vbJrn = mean(abs(dipvbn(v)));
                    mneJrn = mean(abs(dipmnen(v)));
                    Jrn = mean(abs(dipn(v)));

                    Zmsn(r) = Zmsn(r) + msJrn^2;
                    Zvbn(r) = Zvbn(r) + vbJrn^2;
                    Zmnen(r) = Zmnen(r) + mneJrn^2;
                    Zrn(r) = Zrn(r) + Jrn^2;
    %                 if (r == 1)
    %                     Znb = Znb + J^2;
    %                     nnn = nnn + 1;
    %                 end

                    msn(k,r) = msn(k,r) + msJrn * J;
                    vbn(k,r) = vbn(k,r) + vbJrn * J;
                    mnen(k,r) = mnen(k,r) + mneJrn * J;
                    jrn(k,r) = jrn(k,r) + Jrn * J;

                    % error
                    %J = dipn(j);
                    errms(k,r) = errms(k,r) + (J - msJrn)^2;
                    errvb(k,r) = errvb(k,r) + (J - vbJrn)^2;
                    errmne(k,r) = errmne(k,r) + (J - mneJrn)^2;

                    % get next neighbors
                    %v = vizinhos(j,r);
                end
            end

    %        totalndips=size(dip,1);

            Zms = sqrt(Zms);
            Zvb = sqrt(Zvb);
            Zmne = sqrt(Zmne);
            Zr = sqrt(Zr);

            Zmsn = sqrt(Zmsn);
            Zvbn = sqrt(Zvbn);
            Zmnen = sqrt(Zmnen);
            Zrn = sqrt(Zrn);

            for r = 1:R+1
                errms(k,r) = sqrt(errms(k,r));
                errvb(k,r) = sqrt(errvb(k,r));
                errmne(k,r) = sqrt(errmne(k,r));

                ms(k,r) = ms(k,r) / (Z*Zms(r));
                vb(k,r) = vb(k,r) / (Z*Zvb(r));
                mne(k,r) = mne(k,r) / (Z*Zmne(r));
                jr(k,r) = jr(k,r) / (Z*Zr(r));

                msn(k,r) = msn(k,r) / (Zn*Zmsn(r));
                vbn(k,r) = vbn(k,r) / (Zn*Zvbn(r));
                mnen(k,r) = mnen(k,r) / (Zn*Zmnen(r));
                jrn(k,r) = jrn(k,r) / (Zn*Zrn(r));
            end
            %     toc;

            xt = xt ./ totalndips;
            x = x + xt;

        end
        
        k = k+1;

    end
end

x = x ./ k;

S = {vb ms mne jr vbn msn mnen jrn dvb dms dmne x xt errvb errms errmne vizinhos};

% % x = linspace(1,R+1, R+1)
% x = x ./ n;
% 
% A = zeros(1,5);
% for i=1:5
%     A(i) = size(ico(i).wh.faces,1);
% end
% A = wh.areas ./ (100*[A size(wh.orig.faces,1)]);
% 
% x = x * A(o);
% 
% vbm = mean(vb);
% msm = mean(ms);
% mnem = mean(mne);
% jrm = mean(jr);

% figure;
% plot(x, log(msm), '.r');
% hold all;
% plot(x, log(vbm), '.b');
% plot(x, log(mnem), '.g');
% plot(x, log(jrm), '.k');
% 
% legend ('\rho MS', '\rho VB', '\rho MNE', '\rho Real');
% xlim([x(1)-1 x(R+1)+1])
% ylabel '\rho'
% xlabel 'Area cm^2'

% vbnm = mean(vbn);
% msnm = mean(msn);
% mnenm = mean(mnen);
% jrnm = mean(jrn);
% 
% figure;
% plot(x, log(msnm), '.r');
% hold all;
% plot(x, log(vbnm), '.b');
% plot(x, log(mnenm), '.g');
% plot(x, log(jrnm), '.k');
% legend ('\rho MSN', '\rho VBN', '\rho MNEN', '\rho Real');
% xlim([x(1)-1 x(R+1)+1])
% ylabel 'Normalized \rho'
% xlabel 'Area cm^2'
% 
% fprintf ('%d ', x)
% fprintf ('\n')
% fprintf ('%d ', ceil(xt))
% fprintf ('\n')
