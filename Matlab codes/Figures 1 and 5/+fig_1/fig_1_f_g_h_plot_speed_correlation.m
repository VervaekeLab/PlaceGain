

corr_PM_session = []; corr_NM_session = []; corr_ns_session = [];
corr_PM = []; corr_NM = []; corr_ns = [];
allSignal = [];
allSignal_PM = []; allSignal_NM = []; allSignal_ns = [];
for i = 1:length(sessions)
    corr_PM = [corr_PM  speedCorr{i}.val(speedCorr{i}.val > pos_chance)];
    corr_NM = [corr_NM speedCorr{i}.val(speedCorr{i}.val < neg_chance)];
    corr_ns = [corr_ns speedCorr{i}.val((speedCorr{i}.val < pos_chance) & (speedCorr{i}.val > neg_chance))];
       
    corr_PM_session(i) = length(speedCorr{i}.val(speedCorr{i}.val > pos_chance))./length(speedCorr{i}.val);
    corr_NM_session(i) = length(speedCorr{i}.val(speedCorr{i}.val < neg_chance))./length(speedCorr{i}.val);
    corr_ns_session(i) = length(speedCorr{i}.val((speedCorr{i}.val < pos_chance) & (speedCorr{i}.val > neg_chance)))./length(speedCorr{i}.val);;
    
    allSignal_PM{i}  =  mData(i).binnedTrialedSingleROI(speedCorr{i}.val > pos_chance,:);
    allSignal_NM{i}  =  mData(i).binnedTrialedSingleROI(speedCorr{i}.val < neg_chance,:);
    allSignal_ns{i}  =  mData(i).binnedTrialedSingleROI((speedCorr{i}.val < pos_chance) & (speedCorr{i}.val > neg_chance),:);
end


fig = figure();
% subplot(1,2,1)
X = [length(corr_PM), length(corr_NM),length(corr_ns)];
perc = round(100*(X/sum(X)),1);
labels = {strcat(num2str(perc(1)),'%'),strcat(num2str(perc(2)),'%'),strcat(num2str(perc(3)),'%')};
p = pie(X, labels);
set(p(2), 'FontName', 'Gotham', 'FontSize',15)
set(p(4), 'FontName', 'Gotham', 'FontSize',15)
set(p(6), 'FontName', 'Gotham', 'FontSize',15)

set(p(1), 'Facecolor',[72,107, 142]./257, 'EdgeColor', 'none');
set(p(3), 'FaceColor',[68 195 214]./257, 'EdgeColor', 'none');
set(p(5), 'Facecolor',[0.7  0.7 0.7], 'EdgeColor', 'none');

saveas(fig,fullfile(save_dir,'pie_chart'), 'fig')
saveas(fig,fullfile(save_dir,'pie_chart'), 'png')

fprintf('%.2f +/- %.2f %% of cells are positively speed modulated (mean +/- std across sessions) \n',nanmean(100*corr_PM_session),nanstd(100*corr_PM_session))
fprintf('%.2f +/- %.2f %% cells are negatively speed modulated (mean +/- std across sessions) \n',nanmean(100*corr_NM_session),nanstd(100*corr_NM_session))
fprintf('%.2f +/- %.2f %% cells are not speed modulated (mean +/- std across sessions) \n',nanmean(100*corr_ns_session),nanstd(100*corr_ns_session))


fig = figure();
% subplot(1,2,1)
X = [nanmean(100*corr_PM_session),nanmean(100*corr_NM_session),100-(nanmean(100*corr_PM_session)+nanmean(100*corr_NM_session))];
perc = round(100*(X/sum(X)),1);
labels = {strcat(num2str(nanmean(100*corr_PM_session)),'%'),strcat(num2str(nanmean(100*corr_NM_session)),'%'),strcat(num2str(100-(nanmean(100*corr_PM_session)+nanmean(100*corr_NM_session))),'%')};
p = pie(X, labels);
set(p(2), 'FontName', 'Gotham', 'FontSize',15)
set(p(4), 'FontName', 'Gotham', 'FontSize',15)
set(p(6), 'FontName', 'Gotham', 'FontSize',15)

set(p(1), 'Facecolor',[239, 101,96]./255, 'EdgeColor', 'none');
set(p(3), 'FaceColor',[74 142 197]./255, 'EdgeColor', 'none');
set(p(5), 'Facecolor',[194 194 194]./255, 'EdgeColor', 'none');

saveas(fig,fullfile(save_dir,'pie_chart_mean_across_sessions'), 'fig')
saveas(fig,fullfile(save_dir,'pie_chart_mean_across_sessions'), 'pdf')

%% average activity for PM/NM/ns
for i = 1:length(sessions)
    sessionmean_PM(i,:) = normalize(nanmean(allSignal_PM{i}),'range', [0 1]);   
    sessionmean_NM(i,:) = normalize(nanmean(allSignal_NM{i},1),'range', [0 1]); 
    sessionmean_ns(i,:) = normalize(nanmean(allSignal_ns{i}),'range', [0 1]); 
end


fig=figure();
F = 1:78;
amean = nanmean(sessionmean_PM);
astd = nanstd(sessionmean_PM)./sqrt(18);
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[72,107, 142]./257,'linestyle','none','FaceAlpha', 0.5);
hold on
plot(amean,'Color', [72,107, 142]./257, 'LineWidth',1.5)

amean = nanmean(sessionmean_NM);
astd = nanstd(sessionmean_NM)./sqrt(18);
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[68 195 214]./257,'linestyle','none','FaceAlpha', 0.5);
hold on
plot(amean,'Color',[68 195 214]./257, 'LineWidth',1.5)


amean = nanmean(sessionmean_ns);
astd = nanstd(sessionmean_ns)./sqrt(18);
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0.7  0.7 0.7],'linestyle','none','FaceAlpha', 0.5);
hold on
plot(amean,'Color',[0.7  0.7 0.7], 'LineWidth',1.5)

ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

ylim([0 1])
xlabel('Position (cm)')
xticks([0 20 40 60 78])
xticklabels({'0' ,'40','80','120','157'})
ylabel('Normalized mean DF/F')


saveas(fig,fullfile(save_dir,'mean_signals'), 'fig')
saveas(fig,fullfile(save_dir,'mean_signals'), 'png')



fig = figure();
hx = histogram(corr_PM,'Facecolor',[72,107, 142]./257, 'BinWidth', 0.024);
hold on
bin_Width = hx.BinWidth;
histogram(corr_NM, 'Facecolor',[68 195 214]./257, 'BinWidth', bin_Width);
histogram(corr_ns, 'Facecolor',[0.7 0.7  0.7], 'BinWidth', bin_Width);
set(gca, 'FontName', 'Gotham', 'FontSize',15)
xlabel('Peak velocity cross correlation')
ylabel('# Cells')
xticks([-0.8 -0.4 0  0.4 0.8])
box off
scatter(nanmean(corr_PM),40,60,[72,107, 142]./257,'filled')
errorbar(nanmean(corr_PM),40,nanstd(corr_PM),'horizontal','Color',[72,107, 142]./257,'LineWidth',1.5)
scatter(nanmean(corr_NM),40,60,[68 195 214]./257,'filled')
errorbar(nanmean(corr_NM),40,nanstd(corr_NM),'horizontal','Color',[68 195 214]./257,'LineWidth',1.5)

fprintf(' Pearson correlation of positively modulated cells and speed is %.2f +/- %.2f (mean +/- std across sessions) \n',nanmean(corr_PM),nanstd(corr_PM))
fprintf(' Pearson correlation of negatively modulated cells and speed is %.2f +/- %.2f (mean +/- std across sessions) \n',nanmean(corr_NM),nanstd(corr_NM))
