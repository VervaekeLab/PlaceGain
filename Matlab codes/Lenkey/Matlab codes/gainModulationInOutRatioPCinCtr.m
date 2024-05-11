function sData = gainModulationInOutRatioPCinCtr(sData)

% calcuting effect of VIP modulation inside the place field (In) and outside the place field (Out).

sData.gainModulationPCinCtr.InOutMeanActRatio = struct;

nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells;
PFStartBin = NaN(nROIs,1);
PFLength = NaN(nROIs,1);
PFPeakAmpl = NaN(nROIs,1);
PFPeakPos = NaN(nROIs,1);
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
                
% generate a common PlaceFieldStartBin and Length Array for the three protocols (some cells are place cells in one, but not in another protocol. I want to include them)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (which bin, which protocol) , search for peak in position tuning curve
    % Opto-off
    PlaceFieldStartBinOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOff = PlaceFieldStartBinOffPre(~isnan(PlaceFieldStartBinOffPre));
    PlaceFieldLengthOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOff = PlaceFieldLengthOffPre(~isnan(PlaceFieldLengthOffPre));
    Peak = NaN(length(PlaceFieldLengthOff),1);
    for m = 1:1:length(PlaceFieldLengthOff) % if more then on place fields occurs
        Peak(m) = max(PosTuningOff(i,PlaceFieldStartBinOff(m):PlaceFieldStartBinOff(m)+PlaceFieldLengthOff(m)-1));
    end
    PeakOff = max(Peak);
    PeakStartOff = PlaceFieldStartBinOff(find(PeakOff));
    PeakPosOff = min(PosTuningOff(find(PosTuningOff(i,1:80)==PeakOff)));
    % Opto-on
    PlaceFieldStartBinOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOn = PlaceFieldStartBinOnPre(~isnan(PlaceFieldStartBinOnPre));
    PlaceFieldLengthOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOn = PlaceFieldLengthOnPre(~isnan(PlaceFieldLengthOnPre));
    Peak = NaN(length(PlaceFieldStartBinOn),1);
    for m = 1:1:length(PlaceFieldStartBinOn) % if more then on place fields occurs
        Peak(m) = max(PosTuningOn(i,PlaceFieldStartBinOn(m):PlaceFieldStartBinOn(m)+PlaceFieldLengthOn(m)-1));
    end
    PeakOn = max(Peak);
    PeakStartOn = PlaceFieldStartBinOn(find(PeakOn));
    %PeakPosOn = 
    % Opto-after
    PlaceFieldStartBinAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinAfter = PlaceFieldStartBinAfterPre(~isnan(PlaceFieldStartBinAfterPre));
    PlaceFieldLengthAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,:);
    PlaceFieldLengthAfter = PlaceFieldLengthAfterPre(~isnan(PlaceFieldLengthAfterPre));
    Peak = NaN(length(PlaceFieldStartBinAfter),1);
    for m = 1:1:length(PlaceFieldStartBinAfter) % if more then on place fields occurs
        Peak(m) = max(PosTuningAfter(i,PlaceFieldStartBinAfter(m):PlaceFieldStartBinAfter(m)+PlaceFieldLengthAfter(m)-1));
    end
    PeakAfter = max(Peak);
    PeakStartAfter = PlaceFieldStartBinAfter(find(PeakAfter));
    % compare place field peaks
    OffOnAfterPeak = [PeakOn PeakOff PeakAfter];
    if max(OffOnAfterPeak) == PeakOff
        PFStartBin(i) = PeakStartOff;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,find(PlaceFieldStartBinOffPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakOff;
        PFPeakPos(i) = PeakStartOff;
    elseif max(OffOnAfterPeak) == PeakOn
        PFStartBin(i) = PeakStartOn;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,find(PlaceFieldStartBinOnPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakOn;
        PFPeakPos(i) = PeakStartOn;
    elseif max(OffOnAfterPeak) == PeakAfter    
        PFStartBin(i) = PeakStartAfter;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,find(PlaceFieldStartBinAfterPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakAfter;
        PFPeakPos(i) = PeakStartAfter;
    end
end
sData.gainModulationPCinCtr.PFStartBin = PFStartBin;
sData.gainModulationPCinCtr.PFLength = PFLength;
sData.gainModulationPCinCtr.PFPeakPos = PFLength;
sData.gainModulationPCinCtr.note = 'PFStartBin and Length: Place cells in opto-off (ctr) protocol, searched for largest peak among PFs';

% InOutMeanActRatio
for k = 1:1:3
     meanActInPlaceField = NaN(nROIs,1); % mean Ca activity in the PF  
     peakActInPlaceField = NaN(nROIs,1); % peak amplitude in PF (mean)
     meanActOutPlaceField = NaN(nROIs,1); % mean Ca activity outside of the PF
     meanActInOutActRatio = NaN(nROIs,1); % ratio of mean activity inside/outside the PF
     peakActInMeanOutRatio = NaN(nROIs,1); % ratio of peak activity inside/ mean activity outside the PF
     BinPosPeakActInPlaceField = NaN(nROIs,1); % position (bin) of peak activity within PF (mean)
     if k == 1
         PosTuning = PosTuningOff;
     elseif k == 2
         PosTuning = PosTuningOn;
     elseif k == 3 
         PosTuning = PosTuningAfter;
     end
    % repeat for all three protocols
    % calculate mean-in-field/mean-out-of-field activity in opto-off vs opto-on trials
    for j = 1:1:length(PlaceCells)
        i = PlaceCells(j);
        meanActInPlaceField(i) = mean(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));  
        peakActInPlaceField(i) = max(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));
        BinPosPeakActInPlaceField(i) = find(PosTuning(i,:)==peakActInPlaceField(i),1);
        meanActOutPlaceField(i) = mean(PosTuning(i,PFStartBin(i)+PFLength(i):(PFStartBin(i)+nBins)-1));  
        meanActInOutActRatio(i) = meanActInPlaceField(i)/meanActOutPlaceField(i);
        peakActInMeanOutRatio(i) = peakActInPlaceField(i)/meanActOutPlaceField(i);
    end
    if k == 1
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField; 
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 2
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 3
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio = peakActInMeanOutRatio;
    end
end

% generate summary
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActInPlaceField;

sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.peakActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.peakActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.peakActInPlaceField;

sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField;

sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActOutPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActOutPlaceField;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActOutPlaceField;

sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio;
sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio;

sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter = NaN(nROIs,3); 
sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOff.peakInMeanOutRatio;
sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoOn.peakInMeanOutRatio;
sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio;


% means
sData.gainModulationPCinCtr.InOutMeanActRatio.means = struct;
sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,3));

sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));

sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOff = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOn= nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOff = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOn = nanmean(sData.gainModulationPCinCtr.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOnPerInOff = sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOn/sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOff;
sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOnPerOutOff = sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOn/sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOff;

% Save file to same path where other files can be found 
 save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end