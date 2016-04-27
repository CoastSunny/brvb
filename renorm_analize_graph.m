function renorm_analize_graph(S, label, bins)

%     fullscreen = get(0,'ScreenSize');

    n = size(S,2);
    
    R = size(S{1}{1},2) - 1;

    at = 14.8408;
    x = at*S{1}{12};
    
    vbm = mean(S{1}{1});
    msm1 = mean(S{1}{2});
%     msm2 = mean(S{2}{2});
%     msm3 = mean(S{3}{2});
%     msm4 = mean(S{4}{2});
%     mnem = mean(S{1}{3});
    jrm = mean(S{1}{4});

    for j=1:2
        
        if j==1
            mycolor = [0 0 1 ; 0 1 0 ; 1 0 0];
        else
            mycolor = [0.5 0.5 0.5; 0.99 0.99 0.99; 0 0 0 ];
        end
        
        fh = figure('Units','Pixels',...
         'Position', [10 10 300 1170],...
         'PaperPositionMode','Auto');

        fsize = 8;
     
        subplot(3,1,1);
        hold all;
%         plot(x, log(msm1), '-o', 'Color', floor(mycolor(3,:)));
%         plot(x, log(msm2), ':o', 'Color', floor(mycolor(3,:)));
%         plot(x, log(msm3), '--o', 'Color', floor(mycolor(3,:)));
%         plot(x, log(msm4), '-.o', 'Color', floor(mycolor(3,:)));
%         plot(x, log(vbm), '-s', 'Color', floor(mycolor(1,:)));
%         plot(x, log(mnem), '-^', 'Color', floor(mycolor(2,:)));
%         plot(x, log(jrm), '-xk');
% 
%         h_l = legend ('\rho MGVB 1', '\rho MGVB 2', '\rho MGVB 3', '\rho MGVB 4', ...
%             '\rho SGVB', ...
%             '\rho MNE', ...
%             '\rho Real');
%         set(h_l, 'FontSize', fsize-2)
%         lgpos = get(h_l, 'Position');
%         nl = lgpos + [0.05 -0.03 0.05 -0.025];
%         set(h_l, 'Position', nl); 
% 
%         xlim([x(1)-1 x(R+1)+1])
% 
%         ylabel ('{\rho}','FontSize', fsize);
% %         ylabh = get(gca,'YLabel');
% %         set(ylabh,'Position',get(ylabh,'Position') + [0 .2 0]);
% 
%         xlabel ('mm^2','FontSize', fsize-1);
%         title ('\rho for different algorithms','FontSize', fsize);
%         set(gca,'FontSize', fsize);


    i=1;
    
    %         vbm = mean(S{i}{5});
    %         msm = mean(S{i}{6});
    %         mnem = mean(S{i}{7});
    %         jrm = mean(S{i}{8});
    % 
    %         f = figure;
    %         plot(x, log(msm), '-sr');
    %         hold all;
    %         plot(x, log(vbm), '-*b');
    %         plot(x, log(mnem), '-dg');
    %         plot(x, log(jrm), '-ok');
    % 
    %         legend ('\rho MS', ... num2str(errms(1))
    %             '\rho VB', ... num2str(errvb(1))
    %             '\rho MNE', ... num2str(errmne(1)) 
    %             '\rho Real');
    %         xlim([x(1)-1 x(R+1)+1])
    %         ylabel '\rho'
    %         xlabel 'Area cm^2'
    %         title ([label ' ' num2str(i)])
    % 
    %         maxfigsize(f);

    %         errvb = sum(cell2mat(S{i}{14}).^2)/N;
    %         errms = sum(cell2mat(S{i}{15}).^2)/N;
    %         errmne = sum(cell2mat(S{i}{16}).^2)/N;
    % 
    %         errvb = mean(S{i}{14});
    %         errms = mean(S{i}{15});
    %         errmne = mean(S{i}{16});
    % 
    %         fe = figure;
    %         plot(x, errms, '-or');
    %         hold all;
    %         plot(x, errvb, '-sb');
    %         plot(x, errmne, '-^g');
    % 
    %         legend (['\epsilon MG ' num2str(i)], ... num2str(errms(1))
    %             '\epsilon VB', ... num2str(errvb(1))
    %             '\epsilon MNE');
    %         xlim([x(1)-1 x(R+1)+1])
    %         ylabel '\epsilon'
    %         xlabel 'Area s (cm^2)'
    %         title (['Erro absoluto ' label ' ' num2str(i)])
    %         
    %         maxfigsize(fe);

            %hist([S{i}{9}(:,1) S{i}{10}(:,1) S{i}{11}(:,1)], 20);

            
            subplot(3,1,2);
            hist([S{i}{9}(:,1) S{i}{11}(:,1) S{i}{10}(:,1)], bins);
            h_legend = legend ('SGVB', 'MNE', ['MGVB ' num2str(i)]);
            set(h_legend,'FontSize',fsize-2);
            lgpos = get(h_legend, 'Position');
            nl = lgpos + [0.03 -0.015 0 0];
            set(h_legend, 'Position', nl); 
        
            title (['Distance between strongest dipoles : ' label ' ' num2str(i)], 'FontSize',fsize);
            set(gca,'FontSize', fsize);
            xlabel ('mm','FontSize', fsize-1);
%             xlim([-3 104]);
            xlim(bins([1,end])+[-10 +10])
            ylim([0 60]);
            colormap(mycolor);
            

            subplot(3,1,3);
            hist([S{i}{9}(:,2) S{i}{11}(:,2) S{i}{10}(:,2)], bins);
            h_legend = legend ('SGVB', 'MNE', ['MGVB ' num2str(i)]);
            set(h_legend,'FontSize',fsize-2);
            lgpos = get(h_legend, 'Position');
            nl = lgpos + [0.03 -0.015 0 0];
            set(h_legend, 'Position', nl); 

            title (['Distance between second strongest dipoles : ' label ' ' num2str(i)], 'FontSize',fsize);
            set(gca,'FontSize', fsize);

            xlabel ('mm','FontSize', fsize-1);
%             xlim([-3 104]);
            xlim(bins([1,end])+[-10 +10])
            ylim([0 60]);
            colormap(mycolor);

% 
    %         f4 = figure;
    %         hist([S{i}{9}(:,3) S{i}{10}(:,3) S{i}{11}(:,3)]);
    % %         h = findobj(gca,'Type','patch');
    % %         f = get(h, 'faces');
    % %         size(f)
    % %         set(h,'FaceColor','flat', 'FaceVertexCData', [0.1; 0.2; 0.9])
    %         legend ('VB', 'MS', 'MNE');
    %         title (['Média das distancias entre dipolos por ordem de intensidade : ' label ' ' num2str(i)]);
    % 

%         print(fh, '-r300', '-dtiff', ['rhos' num2str(j)]);
%         print(fh, '-r300', '-deps', ['rhos' num2str(j)]);
%         print(fh, '-r300', '-dpdf', ['rhos' num2str(j)]);
%         print(fh, '-r300', '-dpng', ['rhos' num2str(j)]);
        
    end
end
