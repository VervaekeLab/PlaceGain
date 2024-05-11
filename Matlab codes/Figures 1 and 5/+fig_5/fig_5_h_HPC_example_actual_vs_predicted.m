
 %   VIP-Arch example
experiment = sessionInfo(4);
session    = 5;

load(fullfile(dirct,experiment.type,experiment.sessions{session},'d_data_off.mat'));
predicted_positon_off = d_data.iter{d_data.medianTestFold}.predPos_test;
real_position_off = d_data.iter{d_data.medianTestFold}.realPos_test;

% load(fullfile(dirct,experiment.type,experiment.sessions{session},'d_data_opto.mat'));
% predicted_positon_opto = d_data.iter{d_data.medianTestFold}.predPos_test;
% real_position_opto = d_data.iter{d_data.medianTestFold}.realPos_test;

fig        = figure('Units','normalized','Position',[0 0 0.5 0.25],'color','w');

ax_off     =  axes('Units','normalized','Position', [0.05 0.2 0.9 0.7]);
% ax_opto    =  axes('Units','normalized','Position', [0.05 0.55 0.9 0.5]);


plot(ax_off,1:length(real_position_off),real_position_off, 'LineWidth', 2.5,'Color','k' )
hold on
plot(ax_off,1:length(predicted_positon_off),predicted_positon_off , 'LineWidth', 2.5, 'Color','b')

% plot(ax_opto,1:length(real_position_off),real_position_off, 'LineWidth', 2,'Color','k' )
% hold on
% plot(ax_opto,1:length(predicted_positon_off),predicted_positon_off , 'LineWidth', 2, 'Color','o')
yline(min(real_position_off),'--','k','LineWidth',2)
xlim(ax_off,[1 930])
xticks([0:31*10:length(real_position_off)])
xticklabels(0:10:10*length(real_position_off)/(31*10))
xlabel(ax_off,'Time (s)')
ylabel(ax_off,'Decoding error (cm)')
yticks(ax_off,[0 round(max(real_position_off))])
yticklabels(ax_off,[0 157])
set(ax_off, 'FontName', 'Arial');
set(ax_off, 'FontSize' ,18);
box(ax_off,'off')

ax_off.XColor = [0 0 0];
ax_off.YColor = [0 0 0];
ax_off.XLabel.Color = [0 0 0];
ax_off.YLabel.Color = [0 0 0];
ax_off.XAxis.LineWidth = 1.5;
ax_off.YAxis.LineWidth = 1.5;
