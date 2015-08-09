function plot_avgresponseSpectrum(ax, ay, E, dt, T_range, np, unitstr, convf, str, nfolders, rsno, nprof, nvals, nchoose1, nchoose2)

defaultFigureProperties;
profid = {'ATCT','Middle','North','South'};

nax = size(ax,2);
fn = 1./T_range;
G = 9.81;

RSx = NaN(np, 7, nax, nfolders);
RSy = NaN(np, 7, nax, nfolders);
RS_x = NaN(np, 7, nax, nfolders);
RS_y = NaN(np, 7, nax, nfolders);

for j = 1:nfolders
    plotind = find(ismember(nvals,reshape(nchoose1,nprof*3,1))>0);
    for i = plotind
        RSx(:,:,i,j) = CalResponseSpectra_SP(fn, dt{j}, ax{j}(:,i)./9.81, E, G);
        RSy(:,:,i,j) = CalResponseSpectra_SP(fn, dt{j}, ay{j}(:,i)./9.81, E, G);
    end



    plotind2 = find(ismember(nvals,reshape(nchoose2,nprof*3,1))>0);
    for i = plotind2
        RS_x(:,:,i,j) = CalResponseSpectra_SP(fn, dt{j}, ax{j}(:,i)./9.81, E, G);
        RS_y(:,:,i,j) = CalResponseSpectra_SP(fn, dt{j}, ay{j}(:,i)./9.81, E, G);
    end
end

if rsno == 1
%     figure;
%     for i = 1:nprof
%         subplot(2,2,i);
%         plotinds = ismember(nvals,nchoose1(:,i));
%         plot(T_range, convf(1)*squeeze(RSx(:,3,plotinds)),'LineWidth',2); grid on;
%         xlabel('Period [s]'); ylabel(strcat('Pseudo-Spectral Acceleration',' [',unitstr{1},']'));
%         legend('Non-Masing','Masing','Linear','Location','NorthEast');
%         title(strcat(str,', X'));
%     end
%     
%     figure;
%     for i = 1:nprof
%         plotinds = ismember(nvals,nchoose1(:,i));
%         subplot(2,2,i);
%         plot(T_range, convf(1)*squeeze(RSy(:,3,plotinds)),'LineWidth',2); grid on;
%         xlabel('Period [s]'); ylabel(strcat('Pseudo-Spectral Acceleration',' [',unitstr{1},']'));
%         legend('Non-Masing','Masing','Linear','Location','NorthEast');
%         title(strcat(profid{i},': ',str,', Y'));
%     end
% elseif rsno == 2 || rsno == 3 % TEMP - no outcrop
%     figure;
%     for i = 1:nprof
%         subplot(2,2,i);
%         plotinds = ismember(nvals,nchoose2(:,i));
%         
%         plot(T_range, convf(1)*squeeze(RS_x(:,3,plotinds)),'LineWidth',2); grid on;
%         xlabel('Period [s]'); ylabel(strcat('Pseudo-Spectral Acceleration',' [',unitstr{1},']'));
%         legend('Non-Masing','Masing','Linear','Location','NorthEast');
%         title(strcat(profid{i},': ',str,', X'));
%     end
%     
%     figure;
%     for i = 1:nprof
%         plotinds = ismember(nvals,nchoose2(:,i));
%         subplot(2,2,i);
%         plot(T_range, convf(1)*squeeze(RS_y(:,3,plotinds)),'LineWidth',2); grid on;
%         xlabel('Period [s]'); ylabel(strcat('Pseudo-Spectral Acceleration',' [',unitstr{1},']'));
%         legend('Non-Masing','Masing','Linear','Location','NorthEast');
%         title(strcat(profid{i},': ',str,', Y'));
%     end
elseif rsno == 4
    figure;
    for i = 1:nprof
        subplot(2,2,i);
        plotinds = ismember(nvals,nchoose1(:,i));
        plotinds2 = ismember(nvals,nchoose2(:,i));
        plot(T_range, squeeze(mean(RSx(:,3,plotinds,:),4))./squeeze(mean(RS_x(:,3,plotinds2,:),4)),'LineWidth',2); grid on;
        xlabel('Period [s]'); ylabel('Factor');
        legend('Non-Masing','Masing','Linear','Location','NorthEast');
        title(strcat(profid{i},': ',str,', X'));
    end
    
    figure;
    for i = 1:nprof
        subplot(2,2,i);
        plotinds = ismember(nvals,nchoose1(:,i));
        plotinds2 = ismember(nvals,nchoose2(:,i));
        plot(T_range, squeeze(mean(RSy(:,3,plotinds,:),4))./squeeze(mean(RS_y(:,3,plotinds2,:),4)),'LineWidth',2); grid on;
        xlabel('Period [s]'); ylabel('Factor');
        legend('Non-Masing','Masing','Linear','Location','NorthEast');
        title(strcat(profid{i},': ',str,', Y'));
    end
    
end


