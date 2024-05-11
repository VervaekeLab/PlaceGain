
clear all
load(strcat('/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/sData/VIP-GCamp/m8057/m8057-20200625-01/sData.mat'));
% load('/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/sData/VIP-GCAMP/m8057/m8057-20200626-00/sData.mat');

fig             = figure('Units','normalized','Position',[0 0 0.3 0.7],'color','w');
ax_lick         = axes('Units','normalized','Position', [0.05 0.04 0.8 0.1]);
ax_speed        =  axes('Units','normalized','Position', [0.05 0.17 0.8 0.1]);
ax_position     =  axes('Units','normalized','Position', [0.05 0.30 0.8 0.1]);
ax_signal        =  axes('Units','normalized','Position', [0.05 0.43 0.8 0.53]);
 
% VIP -GCamp : m8057-20200625-01 [ 9 17 32 35 23 ]; xlimit = [4500 6500];
cells = [ 9 17 32 35 23 ];%[4     5    11    12    14 10];% [ 9 17 32 35 23 ];%36
% 2 
% Parameters
session = 1;
lineWidth = 1;
dff_color = [0,0,250]/255;%[251,176,59]/255;
ind_start = 16000; % bin in recording
ind_end =21000;%size(sData.imData.roiSignals.dff,2); % bin in recording
signal = sData.imdata.roiSignals(2).dff_LPlight;

c = 1; % init

minVal = []; maxVal = [];
% cells = mData(area,session).TimeB.indices;
% For each of the selected landmark cells
for i = 1:length(cells) %[3,5,4,13]   tcs: [3,9,17], pcs: [11,45]
    cell_index = cells(i);
    c = c+1;

    % Plot activity
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,1,1,60); 
    dffToPlot = normalize(signal(cell_index,:),'range',[0,1]);%  normalize(,'range',[0,1]);%,'gaussian',3);

    dffToPlot = dffToPlot + c;
    
    plot(ax_signal,dffToPlot,'LineWidth',lineWidth,'Color','k')
    hold(ax_signal,'on')
    yticks([])    
    box off
    
    maxVal = max([maxVal max(dffToPlot)]);
    minVal = min([minVal min(dffToPlot)]);
end

xlimit = [4500 6500];
ax_speed.XAxis.Visible = 'off'; ax_position.YAxis.Visible = 'off';
plot(ax_position,sData.behavior.wheelPosDs,'k', 'LineWidth',1);
plot(ax_speed,sData.behavior.runSpeedDs,'k', 'LineWidth',1);

[lickEvents,~,~]  =analyseLickData(sData.daqdata.lickSignal, 92, 15, 15);

lick(1) = sum(lickEvents(1:sData.daqdata.frameIndex(1)));
for i = 1:length(sData.daqdata.frameIndex)-1
    lick(i+1) = sum(lickEvents(sData.daqdata.frameIndex(i)+1:sData.daqdata.frameIndex(i+1)));
end

lickFrequency = ifreq(lick,31);
lickFrequency(lickFrequency > 8) = 8;

plot(ax_lick,lickFrequency,'k', 'LineWidth',1);


set(ax_position, 'FontSize',8, 'FontName', 'Arial','YTickLabels',[],'XTickLabels',[], 'Box', 'off', 'XLim', xlimit);
set(ax_speed, 'FontSize',8, 'FontName', 'Arial','YTickLabels',[],'XTickLabels',[], 'Box', 'off', 'XLim', xlimit);
set(ax_signal, 'FontSize',8, 'FontName', 'Arial','YTickLabels',[],'XTickLabels',[], 'Box', 'off', 'XLim', xlimit, 'YLim', [minVal maxVal]);
set(ax_lick, 'FontSize',8, 'FontName', 'Arial','YTickLabels',[],'XTickLabels',[], 'Box', 'off', 'XLim', xlimit);


axes_position.YLabel.String= 'pos';
axes_speed.YLabel.String= 'speed';
axes_lick.YLabel.String= 'lick';


ax_lick.XTick               = xlimit(1):31*10:xlimit(2);
ax_lick.XTickLabel          = 0:10:10*length(ax_lick.XTick);
ax_lick.XLabel.String       = 'Time (s)';

ax_speed.XAxis.Visible = 'off'; axes_position.YAxis.Visible = 'off';


   
function [lickEvents,lickCount,licErrorPercentage] = analyseLickData(lickSignal,samplingRate,lickLengthThreshold,lickFreqThreshold)
lickStartIndexesRaw = find(diff(lickSignal) == 1)+1;
lickEndIndexesRaw = find(diff(lickSignal) == -1);
% Correct for licks not fully captured at the beginning or end of the recording
if lickSignal(1) == 1
    lickEndIndexesRaw = lickEndIndexesRaw(2:numel(lickEndIndexesRaw));
end
if lickSignal(numel(lickSignal)) == 1
    lickStartIndexesRaw = lickStartIndexesRaw(1:numel(lickStartIndexesRaw)-1);
end
lickLengthsMsRaw = (lickEndIndexesRaw - lickStartIndexesRaw)/(samplingRate/1000); %Convert data to ms
diffLickStart = diff(lickStartIndexesRaw);
lickFreqRaw = samplingRate./diffLickStart;
shortLickIndexes = lickStartIndexesRaw(find(lickLengthsMsRaw < lickLengthThreshold));
tooFastLickIndexes = lickStartIndexesRaw(find(lickFreqRaw > lickFreqThreshold));
lickErrorsIndexes = union(shortLickIndexes,tooFastLickIndexes);
% Filtered lick cueOne indexes
lickStartIndexes = setdiff(lickStartIndexesRaw,lickErrorsIndexes);
lickCount = numel(lickStartIndexes);
lickErrors = numel(lickErrorsIndexes);
licErrorPercentage = lickErrors/numel(lickStartIndexesRaw)*100;
% Create lick event array
lickEvents(1:numel(lickSignal)) = zeros;
lickEvents(lickStartIndexes) = 1;

end

function iFrequency = ifreq(signal,fs)
% returns the instantaneous frequency of a discrete signal (containing 0 and 1) with fs sampling rate
% 
%signal = sData.behavior.signals.lickEvents;
%
iFrequency = nan(size(signal));
if size(signal,1) > size(signal,2)
    signal = signal';
end
events = find(signal == 1);
freques = fs./diff(events);
freques = [freques 0];
iFrequency(signal == 1) = freques;
iFrequency(1) = 0;
iFrequency = fillmissing(iFrequency,'previous'); 
end