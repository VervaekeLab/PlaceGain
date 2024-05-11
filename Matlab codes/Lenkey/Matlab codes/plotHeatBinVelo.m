function [VeloMatrix,MeanVeloBin,TRNuExtended] = SA_plotHeatBinVelo(sData,Ymax) 

%PARAMETERS TO BE SET:
TRNu = sData.behavior.wheelLapImaging;
BinNu = sData.behavior.meta.nBins;
EnterIntoBinExtended = sData.behavior.binning.enterIntoBinIndexExtended;
SampleSpentInBinExtended = sData.behavior.binning.SampleSpentInBinExtendedBins;

% get data
Velocity = sData.behavior.runSpeedDs;

% set trial number for plotting, search last full data trial
for i = 1:1:TRNu
    if isnan(EnterIntoBinExtended(i,size(EnterIntoBinExtended,2))) || EnterIntoBinExtended(i,size(EnterIntoBinExtended,2))==0
        TRNuExtended = i-1; 
        break
    else
        TRNuExtended = TRNu;
    end
end

VeloMatrix = NaN(TRNuExtended,size(EnterIntoBinExtended,2));
%VeloMatrixExtended10 = NaN(TRNuPlot,size(EnterIntoBin,2)); % add extra 10 bins to the end of trials from the next trial
for i = 1:1:TRNuExtended  % rows are trials
    for j = 1:1:size(EnterIntoBinExtended,2)  % columns (distance bins)  
      VeloMatrix(i,j) = mean(Velocity(EnterIntoBinExtended(i,j):(EnterIntoBinExtended(i,j)+SampleSpentInBinExtended(i,j)-1)));
      if isnan(VeloMatrix(i,j)) % if one sample has to be shared by two bins extrapolate value for each bin (otherwise it would be a NaN) 
          VeloMatrix(i,j) = mean(Velocity(EnterIntoBinExtended(i,j):(EnterIntoBinExtended(i,j)+1)));
      end
    end
end
MeanVeloBin = nanmean(VeloMatrix,1);

%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:(TRNuExtended),(VeloMatrix(1:TRNuExtended,1:BinNu))) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
j = colorbar;
colormap(jet);
j.Label.String = 'Velocity (cm/s)';
j.Label.FontSize = 11;
j.TickDirection = 'out'; 
caxis([0 Ymax]); %set limits for color plot, below 1st black, above 2nd white

xlabel('Position on wheel (cm)');
ax = gca;
ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
if TRNuExtended < 50
    yticklabels = 0:5:TRNuExtended;
elseif TRNuExtended >= 50 && TRNuExtended < 200
    yticklabels = 0:10:TRNuExtended;
else
    yticklabels = 0:20:TRNuExtended;
end
yticks = linspace(1, (TRNuExtended), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(sData.sessionInfo.fileID);


end