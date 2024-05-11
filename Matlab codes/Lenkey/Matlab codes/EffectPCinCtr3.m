function sData = EffectPCinCtr3(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects all place cells which are categorized as a place cell in control.
% compare peak amlitude and position of peak in control and opto trials
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols
% option to discard data before bin 10, also discard those place cells
% which has peak in the first 10 bin or categorized as landmark cells

%%% CALCULATE EFFECT SIZE OF THE OPTO STIM ON THE MEAN POS TUNING CURVES
sData.effectOnCellLevel.notes = 'analysis on all place cell which has one peak after bin 10'; 
Folder = 'C:\MATLAB\SAVE';
savePath = fullfile(Folder,sData.sessionInfo.fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% ROI IDs which are PC with one peak:
sData.effectPCinCtrOnePFFinal3.PCwithOnePF = find(~isnan(sData.effectPCinCtrOnePFFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(:,1)));
PlaceCellsPre = sData.effectPCinCtrOnePFFinal3.PCwithOnePF;
DiscardBinsAtBegininng = 10; % discard data in the beginning of track
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;

PFPeakPos = NaN(nROIs,1);
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; 
counter = 0;
for i = 1:1:size(PlaceCellsPre,1) % if more then on place fields occurs
    j = PlaceCellsPre(i);
        PeakAmplOff = max(PosTuningOff(j,:));
        if sum(~isnan(PosTuningOff(j,:)))==0 % if there are only NANs
            continue
        end
        PeakPosOff = min(find(PosTuningOff(j,1:nBins)==PeakAmplOff));
        if PeakPosOff < DiscardBinsAtBegininng % if main peak is in the first bins which is to be discarded, discard this ROI from analysis
            continue
        end
        counter = counter +1;
        PlaceCells(counter,1) = j;
        sData.effectPCinCtrOnePFFinal3.PCwithOnePF_PeakAmplOff(counter,1) = PeakAmplOff;
        sData.effectPCinCtrOnePFFinal3.PCwithOnePF_PeakPosOff(counter,1) = PeakPosOff;
end
sData.effectPCinCtrOnePFFinal3.PCwithOnePFPeakAfter10Bin = PlaceCells;

%%% CALCULATE IF THE EFFECT WAS SIGNIFICANT ON THE CELL LEVEL, USING ALL TRIALS
% use non parametric rank sum test on the peak in the PF in ctr and opto
PlaceCellsSelected = sData.effectPCinCtrOnePFFinal3.PCwithOnePFPeakAfter10Bin;
nROIsSelected = length(PlaceCellsSelected);
sData.effectOnCellLevel.ranksumPConePFAllTrials = struct;
ranksumP = NaN(nROIsSelected,1); % P value for the Unsigned Wilcokon test
ranksumH = NaN(nROIsSelected,1); % H=0 non-sign, H=1 sign indicates that the null hypothesis can be rejected 
OptoOnTrials = sData.imdata.binned.OptoOnTrials;
OptoOffTrials = sData.imdata.binned.OptoOffTrials;
%PeaksInPFOff3bin = NaN(length(sData.behavior.opto.OptoOffTrialsIndices),3);
%PeaksInPFOn3bin = NaN(length(sData.behavior.opto.OptoOnTrialsIndices),3);
CellEffectPeak = NaN(nROIsSelected,1);
nonsignificant = 0;
signIncrease = 0;
signDecrease = 0;
% store data
for i = 1:1:nROIsSelected
    sData.effectOnCellLevel.PeaksInPFOffTrials{i} = NaN(length(OptoOffTrials),1);
    sData.effectOnCellLevel.PeaksInPFOnTrials{i} = NaN(length(OptoOnTrials),1);
end
sData.effectOnCellLevel.PeaksInPFOffMeanROIs = NaN(nROIsSelected,1);
sData.effectOnCellLevel.PeaksInPFOnMeanROIs = NaN(nROIsSelected,1);
sData.effectOnCellLevel.EffectPeakOnPerOff = NaN(nROIsSelected,1);
%%% for calculate percentage of increase and decrease on cell level (no matter if it was significant or not)
Incr = 0; 
Decr = 0;
for i = 1:1:nROIsSelected
    j = PlaceCellsSelected(i);
    PFPeakPos = sData.effectPCinCtrOnePFFinal3.PCwithOnePF_PeakPosOff(i);
    data2 = horzcat(sData.imdata.binned.RoidFF{1, j},sData.imdata.binned.RoidFF{1, j}); %binned Ca activity for each roi
    PeaksInPFOff3bin = data2(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1);
    PeaksInPFOff = mean(PeaksInPFOff3bin,2);
    PeaksInPFOn3bin = data2(OptoOnTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); use 1 or 3 bins?
    PeaksInPFOn = mean(PeaksInPFOn3bin,2);
    %%% calcuate if effect is significant
    [ranksumP(i),ranksumH(i)] = ranksum(PeaksInPFOff,PeaksInPFOn); % h=1 significant, h=0 not sign
    CellEffectPeak(i) = mean(PeaksInPFOn)/mean(PeaksInPFOff); % less than 1 decrease
    if ranksumH(i) == 0
       nonsignificant = nonsignificant + 1; 
    elseif CellEffectPeak(i) > 1 && ranksumH(i) == 1
       signIncrease = signIncrease + 1;
    elseif CellEffectPeak(i) < 1 && ranksumH(i) == 1
       signDecrease = signDecrease + 1;
    end
    % calc if in
    if CellEffectPeak(i) > 1 
        Incr = Incr + 1;
    elseif CellEffectPeak(i) < 1
        Decr = Decr + 1;
    end
    sData.effectOnCellLevel.PeaksInPFOffTrials{i} = PeaksInPFOff;
    sData.effectOnCellLevel.PeaksInPFOnTrials{i} = PeaksInPFOn;
    sData.effectOnCellLevel.PeaksInPFOffMeanROIs(i) = nanmean(PeaksInPFOff);
    sData.effectOnCellLevel.PeaksInPFOnMeanROIs(i) = nanmean(PeaksInPFOn);
    sData.effectOnCellLevel.EffectPeakOnPerOff(i) = nanmean(PeaksInPFOn) / nanmean(PeaksInPFOff);
end

% effect
PercentageCellEffectIncr = Incr / nROIsSelected;
PercentageCellEffectDecr = Decr / nROIsSelected;
sData.effectOnCellLevel.PercCellEffectIncr = PercentageCellEffectIncr;
sData.effectOnCellLevel.PercCellEffectDecr = PercentageCellEffectDecr;

sData.effectOnCellLevel.ranksumPConePFAllTrials.ROIs = PlaceCells;
sData.effectOnCellLevel.ranksumPConePFAllTrials.CellEffectPeakOnPerOff = CellEffectPeak;
sData.effectOnCellLevel.ranksumPConePFAllTrials.ranksumP = ranksumP;
sData.effectOnCellLevel.ranksumPConePFAllTrials.ranksumH = ranksumH;
sData.effectOnCellLevel.ranksumPConePFAllTrials.SignIncrNu = signIncrease;
sData.effectOnCellLevel.ranksumPConePFAllTrials.SignDecrNu = signDecrease;
sData.effectOnCellLevel.ranksumPConePFAllTrials.nonsignificantNu = nonsignificant;
sData.effectOnCellLevel.ranksumPConePFAllTrials.SignIncrPercentage = signIncrease/nROIsSelected*100;
sData.effectOnCellLevel.ranksumPConePFAllTrials.SignDecrPercentage = signDecrease/nROIsSelected*100;
sData.effectOnCellLevel.ranksumPConePFAllTrials.nonsignificantPercentage = nonsignificant/nROIsSelected*100;

% plotting
% colors for pie-chart
newColors = [
    0       0       1   % blue
    1       0       0   % red
    0.7     0.7     0.7];  % grey
% plot percentage of increase/decrease
figure('Color','white')
r = pie([PercentageCellEffectDecr PercentageCellEffectIncr]);
r(1).FaceColor = newColors(1,:); 
r(3).FaceColor = newColors(2,:); 

FileName = strcat(sData.sessionInfo.fileID,'-EffectPieIncrDecr'); 
exportgraphics(gcf,fullfile(savePath,[FileName '.png']),'Resolution',300);
exportgraphics(gcf,fullfile(savePath,[FileName '.pdf']),'Resolution',300);
saveas(gcf,(fullfile(savePath,[FileName '.pdf'])));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot percentage of SIGNIFICANT increase/decrease
figure('Color','white')
p = pie([signDecrease signIncrease nonsignificant]);
p(1).FaceColor = newColors(1,:); 
p(3).FaceColor = newColors(2,:); 
p(5).FaceColor = newColors(3,:); 

FileName = strcat(sData.sessionInfo.fileID,'-EffectPieSIGNIncrDecr'); 
exportgraphics(gcf,fullfile(savePath,[FileName '.png']),'Resolution',300);
exportgraphics(gcf,fullfile(savePath,[FileName '.pdf']),'Resolution',300);
saveas(gcf,(fullfile(savePath,[FileName '.pdf'])));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%%% CALCULATE IF THE EFFECT WAS SIGNIFICANT ON THE CELL LEVEL, USING ONLY
%%% TRIALS WHERE THERE IS ACTIVITY IN PEAK
% maybe good for Chrimson recordings, but not for Arch
% use non parametric rank sum test on the peak in the PF in ctr and opto
PlaceCellsSelected = sData.effectPCinCtrOnePFFinal3.PCwithOnePFPeakAfter10Bin;
nROIsSelected = length(PlaceCellsSelected);
sData.effectOnCellLevel.ranksumPConePFActiveTrials = struct;
ranksumP2 = NaN(nROIsSelected,1); % P value for the Unsigned Wilcokon test
ranksumH2 = NaN(nROIsSelected,1); % H=0 non-sign, H=1 sign indicates that the null hypothesis can be rejected 
OptoOnTrials = sData.imdata.binned.OptoOnTrials;
OptoOffTrials = sData.imdata.binned.OptoOffTrials;
%PeaksInPFOff3bin = NaN(length(sData.behavior.opto.OptoOffTrialsIndices),3);
%PeaksInPFOn3bin = NaN(length(sData.behavior.opto.OptoOnTrialsIndices),3);
CellEffectPeak2 = NaN(nROIsSelected,1);
nonsignificant2 = 0;
signIncrease2 = 0;
signDecrease2 = 0;
for i = 1:1:nROIsSelected
    j = PlaceCellsSelected(i);
    PFPeakPos = sData.effectPCinCtrOnePFFinal3.PCwithOnePF_PeakPosOff(i);
    data2 = horzcat(sData.imdata.binned.RoidFF{1, j},sData.imdata.binned.RoidFF{1, j}); %binned Ca activity for each roi
    PeaksInPFOff3bin = data2(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1);
    PeaksInPFOffPre = mean(PeaksInPFOff3bin,2);
    PeaksInPFOn3bin = data2(OptoOnTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); use 1 or 3 bins?
    PeaksInPFOnPre = mean(PeaksInPFOn3bin,2);
    % check if it was an active trial or not
    meanActRoiOutSidePF = sData.effectPCinCtrOnePFFinal3.meanActOutPFCtrPos(j);
    for k = 1:1:length(PeaksInPFOffPre)
        if PeaksInPFOffPre(k)< meanActRoiOutSidePF
            PeaksInPFOffPre(k) = NaN;
        end
    end
    for k = 1:1:length(PeaksInPFOnPre)
        if PeaksInPFOnPre(k)< meanActRoiOutSidePF
            PeaksInPFOnPre(k) = NaN;
        end
    end
    PeaksInPFOff2 = PeaksInPFOffPre(~isnan(PeaksInPFOffPre));
    PeaksInPFOn2 = PeaksInPFOnPre(~isnan(PeaksInPFOnPre));
    if sum(~isnan(PeaksInPFOn2))==0 % if main peak is in the first bins which is to be discarded, discard this ROI from analysis
         continue
    end
    
    %%% calcuate if effect is significant
    [ranksumP2(i),ranksumH2(i)] = ranksum(PeaksInPFOff2,PeaksInPFOn2); % h=1 significant, h=0 not sign
    CellEffectPeak2(i) = mean(PeaksInPFOn2)/mean(PeaksInPFOff2); % less than 1 decrease
    if ranksumH2(i) == 0
       nonsignificant2 = nonsignificant2 + 1; 
    elseif CellEffectPeak2(i) > 1 && ranksumH2(i) == 1
       signIncrease2 = signIncrease2 + 1;
    elseif CellEffectPeak2(i) < 1 && ranksumH2(i) == 1
       signDecrease2 = signDecrease2 + 1;
    end
end
sData.effectOnCellLevel.ranksumPConePFActiveTrials.ROIs = PlaceCells;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.CellEffectPeakOnPerOff = CellEffectPeak2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.ranksumP = ranksumP2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.ranksumH = ranksumH2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.SignIncrNu = signIncrease2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.SignDecrNu = signDecrease2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.nonsignificantNu = nonsignificant2;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.SignIncrPercentage = signIncrease2/nROIsSelected*100;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.SignDecrPercentage = signDecrease2/nROIsSelected*100;
sData.effectOnCellLevel.ranksumPConePFActiveTrials.nonsignificantPercentage = nonsignificant2/nROIsSelected*100;


%%% CALCULATE IF THE EFFECT WAS SIGNIFICANT ON SESSION LEVEL % added 2023.06.30.
% use non parametric PAIRED rank sum test (Wilcoxon) on the peak in the PF
% in ctr and opto for each cell
PeakOff = sData.effectOnCellLevel.PeaksInPFOffMeanROIs; 
PeakOn = sData.effectOnCellLevel.PeaksInPFOnMeanROIs;
% WPvalue: P value for the Unsigned Wilcokon test
% WH: H=0 non-sign, H=1 sign indicates that the null hypothesis can be rejected 

%CellEffectPeak = NaN(nROIsSelected,1);
[WPvalue,WH] = signrank(PeakOff,PeakOn); 
CellEffectPeak = mean(PeakOn)/mean(PeakOff); % less than 1 decrease

sData.effectOnCellLevel.WilcoxonEffectSession.MeanEffectPeakInPFonCells = CellEffectPeak;
sData.effectOnCellLevel.WilcoxonEffectSession.Pvalue = WPvalue;
sData.effectOnCellLevel.WilcoxonEffectSession.ifSignificant = WH;


%%% CALCULATE ZERO SHIFTED POS TUNING CURVES AND CALCULATE SHIFT AND GAIN FOR EACH CELL
sData.effectOnCellLevel.PosTuningShiftedCells = struct;
BinsAroundPeak = 30;
PosTuningShiftedOff = NaN(nROIsSelected,2*BinsAroundPeak+1);
PosTuningShiftedOn = NaN(nROIsSelected,2*BinsAroundPeak+1);
ShiftEffectSize = NaN(nROIsSelected,1);
GainEffectSize = NaN(nROIsSelected,1);
GainShiftRatio = NaN(nROIsSelected,1);
GainEffectSizePercentage = NaN(nROIsSelected,1);
ShiftEffectSizePercentage = NaN(nROIsSelected,1);

for i = 1:1:nROIsSelected
    j = PlaceCells(i);
    PFPeakPos = sData.effectPCinCtrOnePFFinal3.PCwithOnePF_PeakPosOff(i);
    PTOff = PosTuningOff(j,:);
    PTOn = PosTuningOn(j,:);
    if PFPeakPos > BinsAroundPeak
        PosTuningShiftedOff(i,:) = PTOff((PFPeakPos-BinsAroundPeak):(PFPeakPos+BinsAroundPeak));
        PosTuningShiftedOn(i,:) = PTOn((PFPeakPos-BinsAroundPeak):(PFPeakPos+BinsAroundPeak));
    else
        PosTuningShiftedOff(i,:) = PTOff((PFPeakPos+sData.behavior.meta.nBins-BinsAroundPeak):(PFPeakPos+sData.behavior.meta.nBins+BinsAroundPeak));
        PosTuningShiftedOn(i,:) = PTOn((PFPeakPos+sData.behavior.meta.nBins-BinsAroundPeak):(PFPeakPos+sData.behavior.meta.nBins+BinsAroundPeak));
    end
    % gain moulation and shift calculation for each Place cell ROIs
    AmplChangeinPFPeakOnMinusOff = PosTuningShiftedOn(i,BinsAroundPeak+1) - PosTuningShiftedOff(i,BinsAroundPeak+1);
    AmplChangeBeforePFOnMinusOff = PosTuningShiftedOn(i,BinsAroundPeak-24) - PosTuningShiftedOff(i,BinsAroundPeak-24); %baseline change before the peak
    ShiftEffectSize(i) = abs(AmplChangeBeforePFOnMinusOff);
    GainEffectSize(i) = abs(AmplChangeinPFPeakOnMinusOff - AmplChangeBeforePFOnMinusOff);
    GainShiftRatio(i) = GainEffectSize(i)/ShiftEffectSize(i);
    GainEffectSizePercentage(i) = GainEffectSize(i)/(GainEffectSize(i) + ShiftEffectSize(i));
    ShiftEffectSizePercentage(i) = ShiftEffectSize(i)/(GainEffectSize(i) + ShiftEffectSize(i));
end

sData.effectOnCellLevel.PosTuningShiftedCells.PosTuningShiftedOff = PosTuningShiftedOff;
sData.effectOnCellLevel.PosTuningShiftedCells.PosTuningShiftedOn = PosTuningShiftedOn;
sData.effectOnCellLevel.PosTuningShiftedCells.ShiftEffectSize = ShiftEffectSize;
sData.effectOnCellLevel.PosTuningShiftedCells.GainEffectSize = GainEffectSize;
sData.effectOnCellLevel.PosTuningShiftedCells.GainShiftRatio = GainShiftRatio;
sData.effectOnCellLevel.PosTuningShiftedCells.GainEffectSizePercentage = GainEffectSizePercentage;
sData.effectOnCellLevel.PosTuningShiftedCells.ShiftEffectSizePercentage = ShiftEffectSizePercentage;


% Gain modulation and shift on mean pos tuning curve
% calculate change in peak and in baseline, Peak is in bin 31,  Baseline calculation is at bin 6, 25 bin before peak (50 cm)
sData.effectOnCellLevel.meanPosTuningShifted = struct;

meanPosTuningShiftedOff = nanmean(PosTuningShiftedOff,1);
meanPosTuningShiftedOn = nanmean(PosTuningShiftedOn,1); 
meanAmplChangeinPFPeakOnMinusOff = meanPosTuningShiftedOn(BinsAroundPeak+1) - meanPosTuningShiftedOff(BinsAroundPeak+1);
meanAmplChangeBeforePFOnMinusOff = meanPosTuningShiftedOn(BinsAroundPeak-24) - meanPosTuningShiftedOff(BinsAroundPeak-24); %baseline change before the peak
ShiftEffectSize = abs(meanAmplChangeBeforePFOnMinusOff);
GainEffectSize = abs(meanAmplChangeinPFPeakOnMinusOff - meanAmplChangeBeforePFOnMinusOff);
GainShiftRatio = GainEffectSize/ShiftEffectSize;
GainEffectSizePercentage = GainEffectSize/(GainEffectSize + ShiftEffectSize);
ShiftEffectSizePercentage = ShiftEffectSize/(GainEffectSize + ShiftEffectSize);

sData.effectOnCellLevel.meanPosTuningShifted.PosTuningShiftedOff = PosTuningShiftedOff;
sData.effectOnCellLevel.meanPosTuningShifted.PosTuningShiftedOn = PosTuningShiftedOn;
sData.effectOnCellLevel.meanPosTuningShifted.meanPosTuninShiftedOff = meanPosTuningShiftedOff;    
sData.effectOnCellLevel.meanPosTuningShifted.meanPosTuninShiftedOn = meanPosTuningShiftedOn;   
sData.effectOnCellLevel.meanPosTuningShifted.meanAmplChangeinPFPeakOnvsOff = meanAmplChangeinPFPeakOnMinusOff; 
sData.effectOnCellLevel.meanPosTuningShifted.meanAmplChangeBeforePFPeakOnvsOff = meanAmplChangeBeforePFOnMinusOff;  
sData.effectOnCellLevel.meanPosTuningShifted.ShiftEffectSize = ShiftEffectSize;
sData.effectOnCellLevel.meanPosTuningShifted.GainEffectSize = GainEffectSize;
sData.effectOnCellLevel.meanPosTuningShifted.GainShiftRatio = GainShiftRatio;
sData.effectOnCellLevel.meanPosTuningShifted.GainEffectSizePercentage = GainEffectSizePercentage;
sData.effectOnCellLevel.meanPosTuningShifted.ShiftEffectSizePercentage = ShiftEffectSizePercentage;

%%% calculate Regression lines (linear fit) for two datapoints: 1:(outofPCctr,outofPFopto) 2:(inPFctr,inPFopto)
sData.effectOnCellLevel.RegrForTwoPoints = struct;
EqSteepness = NaN(nROIsSelected,1); % steepness of equation
EqShift = NaN(nROIsSelected,1); % shift of equation
% caclulate 2 point: baseline (ctr,opto) and peak(ctr,opto)
PosTuningOffNaN = PosTuningOff; 
PosTuningOnNaN = PosTuningOn;
PosTuningOffNaN(:,1:DiscardBinsAtBegininng) = NaN; % discard data from first bins
PosTuningOnNaN(:,1:DiscardBinsAtBegininng) = NaN;
PosTuningOffNaN(:,sData.behavior.meta.nBins+1:sData.behavior.meta.nBins+DiscardBinsAtBegininng) = NaN; 
PosTuningOnNaN(:,sData.behavior.meta.nBins+1:sData.behavior.meta.nBins+DiscardBinsAtBegininng) = NaN;
PFStartBin = sData.effectPCinCtrOnePFFinal3.PFStartPosBinOffOnAfter(:,1);
PFLength = sData.effectPCinCtrOnePFFinal3.PFLengthBinOffOnAfter(:,1);
%REGRESSION FOR 2 POINTS: not to overrepresent baseline points
%{
for i = 1:1:nROIsSelected
    roi = PlaceCells(i);
    % activity at PF peak in ctr and opto
    PeakInPFOff3bin = nanmean(PosTuningOffNaN(roi,PFPeakPos-1:PFPeakPos+1)); 
    PeakInPFOn3bin = nanmean(PosTuningOnNaN(roi,PFPeakPos-1:PFPeakPos+1));
    % activity out of PF (baseline
    BaselineBins = NaN(6,1); % data of 6 bins will be averaged
    MinPosTuningOff = min(PosTuningOffNaN(roi,PFStartBin(roi)+ :sData.behavior.meta.nBins+PFStartBin(roi))); % search the minimum in the PosTuning curve expect the place fieald
    PosMinPosTuningOff = min(find(PosTuningOffNaN(roi,:) == MinPosTuningOff)); % position of minimum in opto-off
    BaselineBins(1:3) = [PosMinPosTuningOff-1 PosMinPosTuningOff PosMinPosTuningOff+1];
    MinPosTuningOn = min(PosTuningOnNaN(roi,PFStartBin(roi)+PFLength(roi):sData.behavior.meta.nBins+PFStartBin(roi)));
    PosMinPosTuningOn = min(find(PosTuningOnNaN(roi,:) == MinPosTuningOn)); % position of minimum in opto-on
    BaselineBins(4:6) = [PosMinPosTuningOn-1 PosMinPosTuningOn PosMinPosTuningOn+1];
    BBins = unique(BaselineBins);
    BaselineOutPFOff6Bins = nanmean(PosTuningOffNaN(roi,BBins));
    BaselineOutPFOn6Bins = nanmean(PosTuningOnNaN(roi,BBins));
    
    % calculate linear fit for data for each roi for each bin
    PosTunOffTwoPoints = [BaselineOutPFOff6Bins BaselineOutPFOn6Bins];
    PosTunOnTwoPoints = [PeakInPFOff3bin PeakInPFOn3bin];
    %Max = max(max(PosTunOff),max(PosTunOnTwoPoints));
    [~, ~,~, ~, st, sh] = linFit(PosTunOffTwoPoints,PosTunOnTwoPoints);
    %Regression2(i) = R2;
    %SignificanceP(i) = P;
    EqSteepness(i) = st; % steepness of equation
    EqShift(i) = sh; % shift of equation
end
sData.effectOnCellLevel.RegrForTwoPoints.EqSteepness = EqSteepness; 
sData.effectOnCellLevel.RegrForTwoPoints.EqShift = EqShift; 
sData.effectOnCellLevel.RegrForTwoPoints.meanSteepness = mean(EqSteepness);
sData.effectOnCellLevel.RegrForTwoPoints.meanShift = mean(EqShift);
sData.effectOnCellLevel.RegrForTwoPoints.note = 'all place cell which has peak after bin 10 and not landmark cells, two points were generated for each cell: one ctr-opto pair where the PF peak (3 bin average) and 3+3 bin average of around the minimum in ctr and opto respectively ) is and '; 
%}

%%% calculate Regression lines (linear fit) for all 80  datapoints for each cells
sData.effectOnCellLevel.RegrForAllBins = struct;
EqSteepness = NaN(nROIsSelected,1); % steepness of equation
EqShift = NaN(nROIsSelected,1); % shift of equation
Regression2 = NaN(nROIsSelected,1); % shift of equation
SignificanceP = NaN(nROIsSelected,1); % shift of equation
%REGRESSION FOR ALL POINTS
for i = 1:1:nROIsSelected
    roi = PlaceCells(i);
    % calculate linear fit for data for each roi for each bin
    PosTunOffAllBins = PosTuningOff(roi,DiscardBinsAtBegininng+1:nBins);
    PosTunOnAllBins = PosTuningOn(roi,DiscardBinsAtBegininng+1:nBins);
    [~, ~,r2, p, st, sh] = linFit(PosTunOffAllBins,PosTunOnAllBins);
    Regression2(i) = r2;
    SignificanceP(i) = p;
    EqSteepness(i) = st; % steepness of equation
    EqShift(i) = sh; % shift of equation
end
sData.effectOnCellLevel.RegrForAllBins.EqSteepness = EqSteepness; 
sData.effectOnCellLevel.RegrForAllBins.EqShift = EqShift; 
sData.effectOnCellLevel.RegrForAllBins.Regression2 = Regression2;
sData.effectOnCellLevel.RegrForAllBins.SignificanceP = SignificanceP;
sData.effectOnCellLevel.RegrForAllBins.meanSteepness = mean(EqSteepness);
sData.effectOnCellLevel.RegrForAllBins.meanShift = mean(EqShift);
sData.effectOnCellLevel.RegrForAllBins.meanRegression2 = mean(Regression2);
sData.effectOnCellLevel.RegrForAllBins.meanSignificanceP = mean(SignificanceP);
sData.effectOnCellLevel.RegrForAllBins.note = 'all place cell which has peak after bin 10 and not landmark cells, two points were generated for each cell: one ctr-opto pair where the PF peak (3 bin average) and 3+3 bin average of around the minimum in ctr and opto respectively ) is and '; 


%%% REGRESSION FOR ALL DATAPOINTS IN THE SESSION (PLACE CELL SELECTION)
PosTuningOffOnArray = NaN(nROIsSelected*(sData.behavior.meta.nBins-DiscardBinsAtBegininng),2);
PosTuningOffArrayPre = PosTuningOff(PlaceCells,DiscardBinsAtBegininng+1:sData.behavior.meta.nBins);
PosTuningOnArrayPre = PosTuningOn(PlaceCells,DiscardBinsAtBegininng+1:sData.behavior.meta.nBins);
PosTuningOffOnArray(:,1) = PosTuningOffArrayPre(:);
PosTuningOffOnArray(:,2) = PosTuningOnArrayPre(:);
[~, ~, R2, P, st, sh] = linFit(PosTuningOffOnArray(:,1),PosTuningOffOnArray(:,2));
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.PosTuningOffOnArray = PosTuningOffOnArray;
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.Regression2 = R2;
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.SignificanceP = P;
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.EqSteepness = st; 
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.EqShift = sh; 
sData.effectOnCellLevel.RegrAllDataPointsPooledInSession.note = 'all place cell which has peak after bin 10 and not landmark cells, all individual ctr-opto data pairs from the binned Ca data'; 

% Save file to same path where other files can be found 
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'.mat')),'sData');

end