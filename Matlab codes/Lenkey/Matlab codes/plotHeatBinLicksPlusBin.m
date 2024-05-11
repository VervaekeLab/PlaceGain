function [BinLickMatrixExtended,MeanLickBinExtended,LickCmMatrixExtended,MeanLickCmExtended] =  plotHeatBinLicksPlusBin(sData) 

%PARAMETERS TO BE SET:
TRNu = sData.behavior.wheelLap;
BinNu = sData.behavior.meta.nBins;
EnterIntoBin = sData.behavior.binning.enterIntoBinIndexExtended;
SampleSpentInBin = sData.behavior.binning.SampleSpentInBinExtendedBins;
BinNuPlus = size(EnterIntoBin,2);

% READ DATA FROM FILE
Temp = 0; % insert an extra zero to the beginnning
Licking = vertcat(Temp,sData.behavior.lickDs); % check if there is lick within at least half the duration if the frame scan, if yes, consider lick (=1)

% set trial number for plotting, search last full data trial
for i = 1:1:TRNu
    if any(isnan(EnterIntoBin(i,:))) || any(EnterIntoBin(i,:)==0)
        TRNuPlot = i-1; 
        break
    else
        TRNuPlot = i;
    end
end


BinLickMatrixExtended = NaN(TRNuPlot,BinNuPlus);
% calculate the number of licks during each bin (lick/cm)
for i = 1:1:TRNuPlot  % rows are trials
    for j = 1:1:BinNuPlus  % columns (distance bins)  
        LR = 0; % set lickrate to zero before counting in each bin
        for m = 1:1:SampleSpentInBin(i,j) % Note: if the recorging starts with lick, it will not detect that first lick
            if Licking(EnterIntoBin(i,j) + m - 2) == 0 && Licking(EnterIntoBin(i,j) + m-1) == 1   % search the first timepoint (sample) when licking starts, zero changes to one, this will be a lick event
                LR = LR + 1;  % if there is a lick event, increase lick rate with one
            end
        end
        BinLickMatrixExtended(i,j) = LR; % put the cumulative lick number in a bin into BinLick matrix
    end
end
MeanLickBinExtended = nanmean(BinLickMatrixExtended,1);

LickCmMatrixExtended = BinLickMatrixExtended./sData.behavior.meta.binSize;
MeanLickCmExtended = nanmean(LickCmMatrixExtended,1);
%Cmax = ceil(max(MeanLickCm)*2);

%PLOT FIGURE
figure('Color','white'); 
imagesc(1:sData.behavior.meta.binSize:(BinNu*sData.behavior.meta.binSize),1:TRNuPlot,(LickCmMatrixExtended(1:TRNuPlot,1:BinNu))) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
j = colorbar;
colormap(jet);
j.Label.String = 'Licking (licks/cm)';
j.Label.FontSize = 11;
j.TickDirection = 'out'; 
caxis([0 2]); 
hold on;
line([157 157],[0 TRNuPlot],'Color','white');
hold on
line([163 163],[0 TRNuPlot],'Color','white');

xlabel('Position on wheel (cm)');
ax = gca;
ax.TickDir = 'out';
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
if TRNuPlot < 50
    yticklabels = 0:5:TRNuPlot;
elseif TRNuPlot >= 50 && TRNuPlot < 200
    yticklabels = 0:10:TRNuPlot;
else
    yticklabels = 0:20:TRNuPlot;
end
yticks = linspace(1, TRNuPlot, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(sData.sessionInfo.fileID);

end