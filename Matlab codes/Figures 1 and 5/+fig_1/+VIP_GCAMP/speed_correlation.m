

lag_time_bins = 31;
save_dir = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/analysis/VIP-GCAMP/lag_1s';
speedCorr = cell(length(sessions),1);

% this 
rng(1);
% VIP_GCAMP.chance_speed_correlation();

for m = 1:length(sessions)
    

    signal       = mData(m).sData.imdata.roiSignals(2).dff;
    speed        = smoothdata(normalize(mData(m).sData.behavior.runSpeedDs,'range',[0 1]),'gaussian',10);
    position     = mData(m).sData.behavior.wheelPosDsBinned;
    trials       = mData(m).sData.trials.trialLength;
   

    % make sure all variables have the same length
    minLength = min([length(speed),size(signal,2), length(position)]);
    speed     = speed(1:minLength);
    position  = position(1:minLength);
    signal    = signal(:,1:minLength);
    trials    = trials(1:minLength);
    
    % remove NaN values
    idxDelete = find(isnan(speed));
    signal(:,idxDelete) = [];
    position(idxDelete) = [];
    speed(idxDelete) = [];
    trials(idxDelete) = [];
    
    % create saving struct 
    speedCorr{m}.val = [];
    speedCorr{m}.c = [];
    speedCorr{m}.idx  = [];
    speedCorr{m}.posCorr = [];
    speedCorr{m}.negCorr = [];
    
    for j = 1: size(signal,1)

        smoothSignal                            = smoothdata(signal(j,:),'gaussian',10);%,'gaussian',1);%smoothdata(signal(j,:),'gaussian','SmoothingFactor', 0.85),'range',[0,1]);
        [speedCorr{m}.c(:,j),speedCorr{m}.lags] = xcorr(speed-mean(speed),smoothSignal'-mean(smoothSignal), lag_time_bins,'coeff');
        [~ ,speedCorr{m}.idx(j)]                = max(abs(speedCorr{m}.c(:,j)));

        speedCorr{m}.val(j)                     = speedCorr{m}.c(speedCorr{m}.idx(j),j);
        mData(m).binnedTrialedSingleROI(j,:)    = normalize(nanmean(createHeatMapPos(signal, sData,j,position,trials)),'range',[0 1]);
    end
    

    speedCorr{m}.posCorr = find(speedCorr{m}.val> 0);
    speedCorr{m}.negCorr = find(speedCorr{m}.val < 0);  
   
    savePath =  fullfile(save_dir,sessions{m}(1:5),'/',sessions{m});
    if ~exist(savePath, 'dir'); mkdir(savePath); end
    
    speedCorrSum.speedCorr = speedCorr;
    save(fullfile(savePath,'speedCorrSum'),'-struct', 'speedCorrSum');
    
end



function heatMap = createHeatMapPos(signal, sData,idx,position, trialLength)

minLength = min([length(trialLength) length(position), size(signal,2)]);

trialLength = trialLength(1:minLength);
position    = position(1:minLength);
signal      = signal(:,1:minLength);

heatMap = NaN(numel(sData.trials.trialStart)-1,78);

for i = 1:numel(sData.trials.trialStart)-1
    isInTrial   = trialLength == i;
    trialSignal = signal(idx,isInTrial);
    
    for jj = 1:78
        idxInTrialBin                = position(isInTrial) == jj;
        heatMap(i,jj) = nanmean(trialSignal(idxInTrialBin));
    end
    
end

end