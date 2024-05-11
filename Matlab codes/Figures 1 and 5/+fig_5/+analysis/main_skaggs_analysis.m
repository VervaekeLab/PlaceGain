fig_5.load_session_info()

load_dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/sData/VIP-Chr-HPC/'; %for VIP-Chr and VIP-Arch
save_dirct  = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/skaggs_5cm';
for i = 1:length(sessionInfo)

    experiment = sessionInfo(i);
    
      for j = 1:length(experiment.sessions)
    
        dataDirectory = fullfile(load_dirct,experiment.sessions{j}(1:5),experiment.sessions{j});
        load(fullfile(dataDirectory, 'sData.mat'));
        [optoSignal, afterSignal, offSignal, optoPos, afterPos, offPos,optoLap, afterLap, offLap] = preproces_sData(sData);
      
        off_skaggs   = [];
        opto_skaggs  = [];
        after_skaggs = [];
        for p = 1:size( offSignal,1)
            [activity_map_off, occupancy_map]    = createActivityMap(offLap, offPos,  offSignal(p,:));
            [off_skaggs(p),~]                    = fig_5.analysis.calculate_skaggs(occupancy_map, activity_map_off);
            
            [activity_map_after, occupancy_map]  = createActivityMap(afterLap, afterPos,  afterSignal(p,:));
            [after_skaggs(p),~]                  = fig_5.analysis.calculate_skaggs(occupancy_map, activity_map_after);
            
            [activity_map_opto, occupancy_map]   = createActivityMap(optoLap, optoPos,  optoSignal(p,:));
            [opto_skaggs(p),~]                   = fig_5.analysis.calculate_skaggs(occupancy_map, activity_map_opto);

        end
        
        

        if ~exist(fullfile(save_dirct,experiment.type,experiment.sessions{j}))
           mkdir(fullfile(save_dirct,experiment.type,experiment.sessions{j}));
        end 
       
        save(fullfile(save_dirct,experiment.type,experiment.sessions{j},'off_skaggs.mat'),'off_skaggs')
        save(fullfile(save_dirct,experiment.type,experiment.sessions{j},'after_skaggs.mat'),'after_skaggs')
        save(fullfile(save_dirct,experiment.type,experiment.sessions{j},'opto_skaggs.mat'),'opto_skaggs')

      end
 
end

   
function [optoSignal, afterSignal, offSignal, optoPos, afterPos, offPos,  optoLap, afterLap, offLap] = preproces_sData(sData)
   params.intensity_idx = 1; 

    % Extract signals
    signal = double(sData.imdata.roiSignals(2).dff);

    % Extract position
    position = sData.behavior.wheelPosDs;

    % Extract signalOpto
    signalOpto = sData.daqdata.optoSignal(sData.daqdata.frameIndex);

    % Provide trials
    nTrials = numel(sData.trials.trialStart) - 1;
    if isfield(sData.behavior.opto,'OptoOnTrials')
        optoTrials = sData.behavior.opto.OptoOnTrials;
        afterTrials = sData.behavior.opto.AfterOptoTrials;
        offTrials = sData.behavior.opto.OptoOffTrials;
    else
       optoTrials = sData.behavior.opto.LightOnTrials;
        afterTrials = sData.behavior.opto.AfterLightTrials;
        offTrials = sData.behavior.opto.LightOffTrials;
    end
    % Find trial indices
    offTrialsIdx = find(offTrials == 1);
    afterTrialsIdx = find(afterTrials == 1);
    optoTrialsIdx = find(optoTrials == 1);

    % Make sure that the optogenetic laser was definitely not on during off trials
    for i = 1:length(offTrialsIdx)
        isInTrial = sData.trials.trialLength == offTrialsIdx(i);
        if nanmean(signalOpto(isInTrial)) > 1
            offTrials(offTrialsIdx(i)) = NaN;
        end
    end
    offTrialsIdx = find(offTrials == 1);

    minLength = min([size(signal, 2), length(sData.trials.trialLength)]);
    sData.trials.trialLength = sData.trials.trialLength(1:minLength);

    % Find indices for each experimental state
    offIdx = []; for i = 1:length(offTrialsIdx); offIdx = [offIdx, find(sData.trials.trialLength == offTrialsIdx(i))]; end
    afterIdx = []; for i = 1:length(afterTrialsIdx); afterIdx = [afterIdx, find(sData.trials.trialLength == afterTrialsIdx(i))]; end
    optoIdx = []; for i = 1:length(optoTrialsIdx); optoIdx = [optoIdx, find(sData.trials.trialLength == optoTrialsIdx(i))]; end

    % Extract signals and positions based on indices
    optoSignal = signal(:, optoIdx);
    afterSignal = signal(:, afterIdx);
    offSignal = signal(:, offIdx);
    
    optoSignal(optoSignal<0) = 0;
    afterSignal(afterSignal<0) = 0;
    offSignal(offSignal<0) =0;

    optoPos = position(optoIdx);
    afterPos = position(afterIdx);
    offPos = position(offIdx);
    
    optoLap = sData.trials.trialLength(optoIdx);
    afterLap = sData.trials.trialLength(afterIdx);
    offLap = sData.trials.trialLength(offIdx);

end

function [activity_map, occupancy_map] = createActivityMap(lap, pos, signal)

    params.binSize = 10.4667/2;%(1/1.5);%(1/1.5);%/(1/1.5);%/2;%(1/1.5);%/(1/1.5);
    params.binGaps  = params.binSize/2:params.binSize:157.5; % create evenly sized bins, that consider circularity
    params.deletePosThreshold = 20;
    
    
    signal(:, pos < params.deletePosThreshold) = [];
    lap(:,pos < params.deletePosThreshold) = [];
    pos(pos < params.deletePosThreshold) = [];
    params.binGaps( params.binGaps < (params.deletePosThreshold-params.binSize/2)) = [];
    
    posBinned = zeros(size(pos));
    for t = 1:length(posBinned)
        [~,ind] = min(abs(params.binGaps -pos(t)));
        posBinned(t) = ind;
    end

    unique_laps = unique(lap);
    for i = 1:length(unique_laps)
        isInTrial   = lap == unique_laps(i);
        trialSignal = signal(isInTrial);

        for jj = 1:length(params.binGaps )
            idxInTrialBin           = posBinned(isInTrial) == jj;
            activity_map(i,jj)      = nanmean(trialSignal(idxInTrialBin));
            occupancy_map(i,jj)     = nansum(idxInTrialBin);
        end
    end

end

