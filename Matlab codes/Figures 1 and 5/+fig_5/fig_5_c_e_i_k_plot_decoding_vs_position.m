clear all
fig_5.load_session_info()
dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/decoding/';


offCol  = [0.00,0.45,0.74];
optoCol = [0.85,0.33,0.10];
afterCol = [0 0 0];

classes =[15.7000   26.1668   36.6334   47.1001   57.5669   68.0335   78.5002   88.9669   99.4336  109.9003  120.3670  130.8337  141.3004 151.7671];
for i = 1:length(sessionInfo)
    
    experiment = sessionInfo(i);
    for j = 1:length(experiment.sessions)
        
        
        load(fullfile(dirct,experiment.type,experiment.sessions{j},'d_data_off.mat'));
        
        for r = 1:5
            [posError_offTemp(r,:),std_posError_offTemp(r,:),std_posError_offTemp(r,:),meanError_offTemp(r,:)]= calculatePosError(d_data.iter{r}.realPos_test,d_data.iter{r}.predPos_test);
        end
        posError_off{i}(j,:) = nanmean(posError_offTemp,1); std_posError_off{i}(j,:) = nanmean(std_posError_offTemp,1); meanError_off{i}(j,:) = nanmean(meanError_offTemp,1);
        
        load(fullfile(dirct,experiment.type,experiment.sessions{j},'d_data_opto.mat'));
        for r = 1:5
            [posError_optoTemp(r,:),std_posError_optoTemp(r,:),std_posError_optoTemp(r,:),meanError_optoTemp(r,:)]= calculatePosError(d_data.iter{r}.realPos_test,d_data.iter{r}.predPos_test);
        end
        posError_opto{i}(j,:) = nanmean(posError_optoTemp,1); std_posError_opto{i}(j,:) = nanmean(std_posError_optoTemp,1); meanError_opto{i}(j,:) = nanmean(meanError_optoTemp,1);
        

    end
    

    %% mean across mice
%     fig =figure();
%     set(fig, 'Units', 'centimeters');
%     set(fig, 'Position', [0 0 20 15]);
% 
% %     animal_Chr 
%     for  m = 1:length(unique(experiment.mouse))
%         amean_temp(m,:) = nanmean(meanError_off{i}(experiment.mouse==m,:),1);
%     end
%     
%     F = classes;
%     amean = nanmean(amean_temp,1);
%     astd = nanstd(amean_temp,1)/sqrt(size(amean_temp,1));
%     patch(gca,[F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],offCol,'linestyle','none','FaceAlpha', 0.1);
%     hold(gca, 'on')
%     plot(gca,F,amean, 'LineWidth', 3, 'Color',offCol)
%     
%     for  m = 1:length(unique(experiment.mouse))
%         amean_temp(m,:) = nanmean(meanError_opto{i}(experiment.mouse==m,:),1);
%     end
%     
%     F = classes;
%     amean = nanmean(amean_temp,1);
%     astd = nanstd(amean_temp,1)/sqrt(size(amean_temp,1));
%     patch(gca,[F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],optoCol,'linestyle','none','FaceAlpha', 0.1);
%     hold(gca, 'on')
%     plot(gca,F,amean, 'LineWidth', 3, 'Color',optoCol)
%    
   
%     ax = gca;
%     ax.FontSize = 16;
%     ax.FontName = 'Arial';
%     ax.XColor = [0 0 0];
%     ax.YColor = [0 0 0];
%     ax.XLabel.Color = [0 0 0];
%     ax.YLabel.Color = [0 0 0];
%     xlabel('Position (cm)')
%     xticks([0 50 100 151])
%     xlim([0 152])
%     xticklabels([0 50 100 157])
%     ylabel('Decoding error(cm)')
%     ylim([0 37])

    
    fig = figure();
    F = classes;
    amean = nanmean(meanError_off{i},1);
    astd = nanstd(meanError_off{i},1)/sqrt(size(meanError_off{i},1));
    patch(gca,[F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],offCol,'linestyle','none','FaceAlpha', 0.1);
    hold(gca, 'on')
    plot(gca,F,amean, 'LineWidth', 3, 'Color',offCol)
    
    amean = nanmean(meanError_opto{i},1);
    astd = nanstd(meanError_opto{i},1)/sqrt(size(meanError_opto{i},1));
    patch(gca,[F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],optoCol,'linestyle','none','FaceAlpha', 0.1);
    hold(gca, 'on')
    plot(gca,F,amean, 'LineWidth', 3, 'Color',optoCol)
    
%     amean = nanmean(meanError_after{l},1);
%     astd = nanstd(meanError_after{l},1)/sqrt(size(meanError_after{l},1));
%     patch(gca,[F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],afterCol,'linestyle','none','FaceAlpha', 0.1);
%     hold(gca, 'on')
%     plot(gca,F,amean, 'LineWidth', 3, 'Color',afterCol)
%     
    ax = gca;
    ax.FontSize = 16;
    ax.FontName = 'Arial';
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.XLabel.Color = [0 0 0];
    ax.YLabel.Color = [0 0 0];
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    xlabel('Position (cm)')
    xticks([0 50 100 151])
    xlim([0 152])
    xticklabels([0 50 100 157])
    ylabel('Decoding error(cm)')
    ylim([0 40])
    
end



function [predPos,predPos_error,std_predError,mean_predError] = calculatePosError(actualClass,predClass)
classes =[15.7000   26.1668   36.6334   47.1001   57.5669   68.0335   78.5002   88.9669   99.4336  109.9003  120.3670  130.8337  141.3004 151.7671];

maxval= max(classes);
for i = 1:length(predClass)
    predError(i) = abs(actualClass(i)-predClass(i));
    if predError(i) > maxval/2 && actualClass(i) > predClass(i)
        predError(i) = abs(maxval-actualClass(i)+predClass(i));
    elseif predError(i) > maxval/2 && actualClass(i) < predClass(i)
        predError(i) = abs(maxval-predClass(i)+actualClass(i));
    end
end


predPos = zeros(length(classes),1);
predPos_error = zeros(length(classes),1);
for i = 1: length(classes)
    idx = find(round(actualClass) == round(classes(i)));
    predClassIdx = predClass(idx);
    if classes(i) < mean(classes)
        idxChange1 = find(predClassIdx-classes(i)> mean(classes));
        predClassIdx(idxChange1) = max(classes)-predClassIdx(idxChange1);
    elseif classes(i) > mean(classes)
        idxChange2 = find(predClassIdx-classes(i) < -mean(classes));
        predClassIdx(idxChange2) = classes(i)+predClassIdx(idxChange2);
    end
    predPos(i) = nanmean(predClassIdx);
    predPos_error(i)   = nanstd(predClassIdx);

    mean_predError(i) = nanmean(predError(idx));
    std_predError(i)= nanstd(predError(idx))/sqrt(length(idx));
    
end


end


function [rmse,predError] = plotDecodingPos(predPos, pos, keepClose, directory, filename, maxval)

fig = figure();
subplot(2,1,1);
plot(1:length(pos),pos)
hold on
plot(1:length(predPos),predPos, 'LineWidth', 1, 'Color', [198, 82, 77]./255)
legend('real position','predicted position', 'eastoutside');

xlim([1 length(predPos)])
xlabel('Pos(Samples)')
ylabel('Position in [cm]')
yticks([0 40 80 120 160 ])
yticklabels({'0','40','80','120','160'})
set(gca, 'FontName', 'Gotham');
set(gca, 'FontSize' ,18);

predError = zeros(length(pos),1);
for i = 1:length(pos)
    predError(i) = abs(pos(i)-predPos(i));
    if predError(i) > maxval/2 && pos(i) > predPos(i)
        predError(i) = abs(maxval-pos(i)+predPos(i));
    elseif predError(i) > maxval/2 && pos(i) < predPos(i)
        predError(i) = abs(maxval-predPos(i)+pos(i));
    end
end

rmse = mean(predError);

subplot(2,1,2);
plot(1:length(pos),predError, 'LineWidth',1, 'Color', [49, 54, 83]./255)

xlim([1 length(predPos)])
xlabel('Pos(Samples)')
ylabel('Error in [cm]')
yticks([0 40 80])
yticklabels({'0','40','80'})
set(gca, 'FontName', 'Gotham');
set(gca, 'FontSize' ,18);

% if strcmp(keepClose,'close') == 1
%     saveas(fig,fullfile(directory,filename))
%     close
% elseif strcmp(keepClose, 'closeNotSave') == 1
%     close
% end

end
