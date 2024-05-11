
speed_limit             = 10; 
min_interval_length_slow_speed = 31;

speedModData = cell(length(sessions));
for m = 1:length(sessions)
    smoothedPlaceCellSignal = [];
    signal      = mData(m).sData.imdata.roiSignals(2).dff;
    speed       = mData(m).sData.behavior.runSpeedDs;
    waterReward = mData(m).sData.behavior.waterRewardDs;
    licking     = mData(m).sData.behavior.lickDs;
    
    % delete ROIs that don't look "healthy" (followed Nora's instructions)
    signal(deletedROIs{m},:)   = [];
    speedModData{m}.speed_in_start_interval = [];
    speedModData{m}.speed_in_stop_interval  = [];
    
    % make sure all variables have the same length
    minLength   = min([size(signal,2),length(speed),length(waterReward)]);
    signal      = signal(:,1:minLength);
    speed       = speed(1:minLength);
    waterRewart = waterReward(1:minLength);
    position    = mData(m).sData.behavior.wheelPosDs(2:minLength+1); 
    
    
   %% find slow intervals & start and stop events
   % Search for the interval during which the mouse speed is slower than the 
   % specified speed_limit. Ensure that this interval lasts for at least 
   % min_interval_length_slow_speed. Save the 2-second period before entering 
   % this interval and the 2-second period after as a speed stop event. 
   % Also, consider the 2 seconds before and after leaving this interval as a start event.
    
    speed_in_stop_interval = [];
    speed_in_start_interval = [];
    start_interval = [];
    stop_interval = [];
    slow_interval = [];
    
    l = 1;
    if length(find(speed < speed_limit))> 1
        slowIdx     = find(speed < speed_limit);
        diffIdxVals = diff(slowIdx);
        diffIdx     = [1;find(diff(slowIdx) > min_interval_length_slow_speed)+1];
        
        if length(diffIdx) <  1
            slow_interval{1} = slowIdx(1):slowIdx(end);
            
            if length(slow_interval{1}) > 62
                stop_interval{1}  = slow_interval{1}(1)-62:slow_interval{1}(1)+62;
                start_interval{1} = slow_interval{1}(end)-62:slow_interval{1}(end)+62;
                
                if isempty(find(stop_interval{1}< 0, 1))
                    speedModData{m}.speed_in_stop_interval(1,:)  = speed(stop_interval{1});
                    speedModData{m}.speed_in_start_interval(1,:) = speed(start_interval{1});
                    all_still_Interval{1} = slow_interval{1};  
                end
                
                l = 2;
            end
        else
            for f = 1:length(diffIdx)-1
                slow_interval{f} = slowIdx(diffIdx(f)):slowIdx(diffIdx(f+1)-1);
                
                if length(slow_interval{f}) > 62
                    stop_interval{l}  = slow_interval{f}(1)-62:slow_interval{f}(1)+62;
                    start_interval{l} = slow_interval{f}(end)-62:slow_interval{f}(end)+62;
                    
                    if isempty(find(stop_interval{l} < 0, 1))
                        speedModData{m}.speed_in_stop_interval(l,:)  = speed(stop_interval{l});
                        speedModData{m}.speed_in_start_interval(l,:) = speed(start_interval{l}); 
                        all_still_Interval{l} = slow_interval{f}; 
                        l = l+1;
                    end
                end
            end
        end
    end
    

    %% find slow intervals & start and stop events

    signal_in_stop_interval = [];
    signal_in_start_interval = [];
    single_cell_signal_in_stop_interval = [];
    single_cell_signal_in_start_interval = [];
    restingPool = [];
    runningPool = [];
    runningModulation = [];
    modulation = [];
    start_modulation = [];
    stop_modulation = [];
    numBootstrap = 300;

    %% Read out cell activity in start and stop events
    %  Make sure that there is no water delivery in the stop event
    %  interval. (Since the stop event is always before start event, we don't discard
    %  the stop intervals in which in the corresponding start event
    %  interval there was a water delivery)
    
    for j = 1: size(signal,1)
        
        smoothSignal       =  normalize(smoothdata(signal(j,:),'gaussian',5),'range',[0,1]);
        all_resting_indices = [];
        all_resting_indices_notReward = [];
        isRew = [];
        single_cell_signal_in_stop_interval = []; single_cell_signal_in_start_interval = [];
        for f = 1:l-1
            
            indices_stop    = waterReward(stop_interval{f});
            indices_start   = waterReward(start_interval{f});
            indices_resting = waterReward(all_still_Interval{f});
            
            
            if isempty(find(indices_stop == 1)) && isempty(find(indices_start == 1))
                isNotRew(f) = 1;
                single_cell_signal_in_start_interval = [single_cell_signal_in_start_interval; smoothSignal(start_interval{f})];
                single_cell_signal_in_stop_interval  = [single_cell_signal_in_stop_interval;  smoothSignal(stop_interval{f})];
            elseif isempty(find(indices_stop == 1)) && ~isempty(find(indices_start == 1))
                isNotRew(f) = 1;
                single_cell_signal_in_stop_interval = [single_cell_signal_in_stop_interval; smoothSignal(stop_interval{f})];
            else
                isRew(f) = 1;
            end
            
            if isempty(find(indices_resting == 1))
                all_resting_indices = [all_resting_indices, all_still_Interval{f}];
            end
            
        end
        
        %% calculate significance of stop modulation with bootstrapping
        if ~isempty(single_cell_signal_in_stop_interval)
            signal_in_stop_interval(j,:)    = normalize(nanmean(single_cell_signal_in_stop_interval,1),'range',[0,1]);% normalize(,'range',[0,1]); %nanmean(normalize(rewSingleSignal,'range',[0,1]),1);
            stop_prePostDiff(j)             = nanmean(signal_in_stop_interval(j,62:93), 'all')-nanmean(signal_in_stop_interval(j,31:62),'all');

            allStop = reshape(cell2mat(stop_interval)', 125,length(stop_interval));
            allPreStop = allStop(31:62,:); allPreStop = allPreStop(:);
            allPostStop = allStop(62:93,:); allPostStop = allPostStop(:);
            
            [~,bootsam_preStop]  = bootstrp(numBootstrap,@mean,allPreStop); %all_resting_indices_notReward);  
            [~,bootsam_postStop]  = bootstrp(numBootstrap,@mean,allPostStop);
    
            boot_start_modulation = NaN(1,numBootstrap);
            for i = 1:numBootstrap   
                bootPreStop = allPreStop(bootsam_preStop(:,i));%all_resting_indices_notReward(bootsam_rest(:,i));
                bootPostStop = allPostStop(bootsam_postStop(:,i));
                boot_stop_modulation(i) = nanmean(smoothSignal(bootPostStop))- nanmean(smoothSignal(bootPreStop));
            end
            
            pd = fitdist(boot_stop_modulation','Normal');
            ci = paramci(pd,'Alpha',.01);  
            
            if ci(1,1) < 0 && ci(2,1) < 0
                stop_modulation(j) = -1;
            elseif ci(1,1)<= 0 && ci(2,1) >= 0
                stop_modulation(j) = 0;
            elseif ci(1,1) > 0 && ci(2,1) > 0
                stop_modulation(j) = 1;
            end
        end
        
        %% calculate significance of start modulation with bootstrapping
        if ~isempty(single_cell_signal_in_start_interval)
            if m == 2
                if j == 31
                    a =1;
                end
            end

            signal_in_start_interval(j,:)     = normalize(nanmean(single_cell_signal_in_start_interval,1),'range',[0,1]);% smoothdata(,'gaussian',10)
            start_prePostDiff(j)              = nanmean(signal_in_start_interval(j,62:93), 'all')-nanmean(signal_in_start_interval(j,31:62), 'all');

            allStart = reshape(cell2mat(start_interval)', 125,length(start_interval));
            allPreStart = allStart(31:62,:); allPreStart = allPreStart(:);
            allPostStart = allStart(62:93,:); allPostStart = allPostStart(:);
                        
            [~,bootsam_preStart]   = bootstrp(numBootstrap,@mean,allPreStart); %all_resting_indices_notReward);  
            [~,bootsam_postStart]  = bootstrp(numBootstrap,@mean,allPostStart);
            
            boot_start_modulation = NaN(1,numBootstrap);
            for i = 1:numBootstrap
                bootPreStart  = allPreStart(bootsam_preStart(:,i));%all_resting_indices_notReward(bootsam_rest(:,i));
                bootPostStart = allPostStart(bootsam_postStart(:,i));
                boot_start_modulation(i) = nanmean(smoothSignal(bootPostStart))- nanmean(smoothSignal(bootPreStart));
            end
            
            
            pd = fitdist(boot_start_modulation','Normal');
            ci = paramci(pd,'Alpha',.01);    

            % A cell was considered to have a significant positive running modulation if 
            % CI is greater than 0, negative modulation if CI is less than 0, and no modulation if the CI contains 0.

            if ci(1,1) < 0 && ci(2,1) < 0
                start_modulation(j) = -1;
            elseif ci(1,1)<= 0 && ci(2,1) >= 0
                start_modulation(j) = 0;
            elseif ci(1,1) > 0 && ci(2,1) > 0
                start_modulation(j) = 1;
            end
        end
        
        % A cell was considered to have a significant positive running modulation if
        % CI is greater than 0, negative modulation if CI is less than 0, and no modulation if the CI contains 0.
        
         if ~isempty(all_resting_indices)
   
            %% calculataed runningModulation
            %% Bootstrapping to test significance of modulation (PVM or NVM)
            % should be in speed_lag_correlation.m, ... but isn't
            restingPool(j,:)     = smoothSignal(all_resting_indices);%_notReward);
            all_runningIndices   = setdiff(1:length(speed), all_resting_indices);
            runningPool(j,:)     = smoothSignal(all_runningIndices);
            runningModulation(j) = nanmean(runningPool(j,:))-nanmean(restingPool(j,:));

            [~,bootsam_rest]  = bootstrp(numBootstrap,@mean,all_resting_indices); %all_resting_indices_notReward);
            [~,bootsam_run]  = bootstrp(numBootstrap,@mean,all_runningIndices);
            for i = 1:numBootstrap
                bootRest = all_resting_indices(bootsam_rest(:,i));%all_resting_indices_notReward(bootsam_rest(:,i));
                bootRun = all_runningIndices(bootsam_run(:,i));
                boot_running_modulation(i) = nanmean(smoothSignal(bootRun))- nanmean(smoothSignal(bootRest));
            end

            pd = fitdist(boot_running_modulation','Normal');
            ci = paramci(pd,'Alpha',.01);

            % A cell was considered to have a significant positive running modulation if
            % CI is greater than 0, negative modulation if CI is less than 0, and no modulation if the CI contains 0.

            if ci(1,1) < 0 && ci(2,1) < 0
                modulation(j) = -1;
            elseif ci(1,1)<= 0 && ci(2,1) >= 0
                modulation(j) = 0;
            elseif ci(1,1) > 0 && ci(2,1) > 0
                modulation(j) = 1;
            end
            
         end
        end
        

            %% NOW: Bootstrapping to check significance for stop and start modulation
           
    speedModData{m}.pvm_start_modulation   = [];
    speedModData{m}.pvm_stop_modulation    = [];
    speedModData{m}.nvm_start_modulation   = [];
    speedModData{m}.nvm_stop_modulation    = [];
    
    speedModData{m}.pvm_signal_in_start_interval   = [];
    speedModData{m}.pvm_signal_in_stop_interval    = [];
    speedModData{m}.nvm_signal_in_start_interval   = [];
    speedModData{m}.nvm_signal_in_stop_interval    = [];
    
    speedModData{m}.pvm_reward_modulation = [];
    speedModData{m}.nvm_reward_modulation = [];
    speedModData{m}.pvm_rewSignal = [];
    speedModData{m}.nvm_rewSignal= [];
   
    speedModData{m}.pvm_run_mod_val = [];
    speedModData{m}.nvm_run_mod_val= [];
    speedModData{m}.pvm_start_prePostDiff = [];
    speedModData{m}.pvm_stop_prePostDiff = [];
    speedModData{m}.nvm_start_prePostDiff = [];
    speedModData{m}.nvm_stop_prePostDiff = [];
    
    speedModData{m}.pvm_reward_modulation_val =[];
    speedModData{m}.nvm_reward_modulation_val = [];
      
    
    
    % load correlation values
    load(strcat(save_dir,'/',sessions{m}(1:5),'/',sessions{m},'/speedCorrSum.mat'))

    speedCorr{m}.val(deletedROIs{m}) = [];
    speedCorr{m}.c(:,deletedROIs{m}) = [];
    
    posCorr = find(speedCorr{m}.val > pos_chance);
    negCorr = find(speedCorr{m}.val < neg_chance);
    
    
    if ~isempty(single_cell_signal_in_start_interval)
        fig = figure();
        subplot(3,4,1)
        plot(nanmean(speedModData{m}.speed_in_start_interval,1),'k', 'LineWidth', 1.5)
        xticks([0 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('velocity (cm/s)')
        xlim([0 124])
        title('speed start')
        hold on
        
        pvm_signal_in_start_interval      = signal_in_start_interval(posCorr,:);
        nvm_signal_in_start_interval      = signal_in_start_interval(negCorr,:);
        pvm_start_modulation = start_modulation(posCorr);
        nvm_start_modulation = start_modulation(negCorr);
        
        pvm_start_prePostDiff = start_prePostDiff(posCorr); [pvm_start_mod_sort_val pvm_start_mod_sort] = sort(pvm_start_prePostDiff);
        nvm_start_prePostDiff = start_prePostDiff(negCorr); [nvm_start_mod_sort_val nvm_start_mod_sort] = sort(nvm_start_prePostDiff);
        speedModData{m}.pvm_signal_in_start_interval = pvm_signal_in_start_interval(pvm_start_mod_sort,:);  speedModData{m}.pvm_start_modulation = pvm_start_modulation(pvm_start_mod_sort);
        speedModData{m}.nvm_signal_in_start_interval = nvm_signal_in_start_interval(nvm_start_mod_sort,:);  speedModData{m}.nvm_start_modulation = nvm_start_modulation(nvm_start_mod_sort);
        
        speedModData{m}.pvmStartCombined = [speedModData{m}.pvm_signal_in_start_interval(speedModData{m}.pvm_start_modulation == 1,:);...
            speedModData{m}.pvm_signal_in_start_interval(speedModData{m}.pvm_start_modulation == 0,:);...
            speedModData{m}.pvm_signal_in_start_interval(speedModData{m}.pvm_start_modulation == -1,:)];
        imAlpha  = ones(size(speedModData{m}.pvmStartCombined));
        imAlpha(isnan(speedModData{m}.pvmStartCombined)) = 0;
        
        subplot(3,4,5)
        imagesc(speedModData{m}.pvmStartCombined, 'AlphaData', imAlpha)
        xticks([1 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('cell')
        title('PVM - start modulation')
        hold on
        xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
        colormap(jet)
        yline(length(find(speedModData{m}.pvm_start_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        yline(length(find(speedModData{m}.pvm_start_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        %
        %     %    caxis([0 1])
        hold on
        set(gca,'ydir','normal');
        
        speedModData{m}.nvmStartCombined = [speedModData{m}.nvm_signal_in_start_interval(speedModData{m}.nvm_start_modulation == 1,:);...
            speedModData{m}.nvm_signal_in_start_interval(speedModData{m}.nvm_start_modulation == 0,:);...
            speedModData{m}.nvm_signal_in_start_interval(speedModData{m}.nvm_start_modulation == -1,:)];
        imAlpha  = ones(size(speedModData{m}.nvmStartCombined));
        imAlpha(isnan(speedModData{m}.nvmStartCombined)) = 0;
        
        subplot(3,4,9)
        imagesc(speedModData{m}.nvmStartCombined, 'AlphaData', imAlpha)
        xticks([1 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('cell')
        title('NVM - start modulation')
        hold on
        xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
        colormap(jet)
        set(gca, 'ydir', 'normal')
        yline(length(find(speedModData{m}.nvm_start_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        yline(length(find(speedModData{m}.nvm_start_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')

        speedModData{m}.pvm_start_prePostDiff = start_prePostDiff(posCorr);
        speedModData{m}.nvm_start_prePostDiff = start_prePostDiff(negCorr);
        
        subplot(3,4,3)
        hobj = bar(1,nanmean(speedModData{m}.pvm_start_prePostDiff), 'FaceColor', [72,107, 142]./257);
        hold on
        errorbar(1, nanmean(speedModData{m}.pvm_start_prePostDiff), nanstd(pvm_start_prePostDiff)/length(pvm_start_prePostDiff), 'k', 'linestyle', 'none', 'LineWidt', 2);
        hobj = bar(2, nanmean(speedModData{m}.nvm_start_prePostDiff), 'FaceColor', [212, 147, 145]./257);
        errorbar(2, nanmean(speedModData{m}.nvm_start_prePostDiff), nanstd(nvm_start_prePostDiff)/length(nvm_start_prePostDiff), 'k', 'linestyle', 'none', 'LineWidt', 2);
        xticks([1 2])
        xticklabels({'PVM', 'NVM'})
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        ylabel('DFF_{post}-DFF_{pre}')
        title('running start')
        ax = gca;
        ax.YAxis.Exponent = 0;
        hold on
    
    end
    
    if ~isempty(single_cell_signal_in_stop_interval)
        subplot(3,4,2)
        plot(nanmean(speedModData{m}.speed_in_stop_interval,1),'k', 'LineWidth', 1.5)
        xticks([0 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('velocity (cm/s)')
        xlim([0 124])
        title('speed stop')
        corrVals = [];
        
        
        pvm_signal_in_stop_interval  = signal_in_stop_interval(posCorr,:);
        nvm_signal_in_stop_interval  = signal_in_stop_interval(negCorr,:);
        
        pvm_stop_modulation  = stop_modulation(posCorr);
        nvm_stop_modulation  = stop_modulation(negCorr);
        
        
        pvm_stop_prePostDiff = stop_prePostDiff(posCorr); [pvm_stop_mod_sort_val pvm_stop_mod_sort] = sort(pvm_stop_prePostDiff);
        nvm_stop_prePostDiff = stop_prePostDiff(negCorr); [nvm_stop_mod_sort_val nvm_stop_mod_sort] = sort(nvm_stop_prePostDiff);
        speedModData{m}.pvm_signal_in_stop_interval = pvm_signal_in_stop_interval(pvm_stop_mod_sort,:);  speedModData{m}.pvm_stop_modulation = pvm_stop_modulation(pvm_stop_mod_sort);
        speedModData{m}.nvm_signal_in_stop_interval = nvm_signal_in_stop_interval(nvm_stop_mod_sort,:);  speedModData{m}.nvm_stop_modulation = nvm_stop_modulation(nvm_stop_mod_sort);
        
        speedModData{m}.pvm_StopCombined = [speedModData{m}.pvm_signal_in_stop_interval(speedModData{m}.pvm_stop_modulation == 1,:);...
            speedModData{m}.pvm_signal_in_stop_interval(speedModData{m}.pvm_stop_modulation == 0,:);...
            speedModData{m}.pvm_signal_in_stop_interval(speedModData{m}.pvm_stop_modulation == -1,:)];
        imAlpha  = ones(size(speedModData{m}.pvm_StopCombined));
        imAlpha(isnan(speedModData{m}.pvm_StopCombined)) = 0;
        
        subplot(3,4,6)
        imagesc(speedModData{m}.pvm_StopCombined, 'AlphaData', imAlpha)
        xticks([1 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('cell')
        title('PVM - stop modulation')
        hold on
        xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
        colormap(jet)
        yline(length(find(speedModData{m}.pvm_stop_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        yline(length(find(speedModData{m}.pvm_stop_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        
        speedModData{m}.nvmStopCombined = [speedModData{m}.nvm_signal_in_stop_interval(speedModData{m}.nvm_stop_modulation == 1,:); ...
            speedModData{m}.nvm_signal_in_stop_interval(speedModData{m}.nvm_stop_modulation == 0,:);...
            speedModData{m}.nvm_signal_in_stop_interval(speedModData{m}.nvm_stop_modulation == -1,:)];
        imAlpha  = ones(size(speedModData{m}.nvmStopCombined));
        imAlpha(isnan(speedModData{m}.nvmStopCombined)) = 0;
        
        subplot(3,4,10)
        imagesc(speedModData{m}.nvmStopCombined,'AlphaData', imAlpha)
        xticks([1 31 62 93 124])
        xticklabels([-2 -1 0 1 2])
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        xlabel('time(s)')
        ylabel('cell')
        title('NVM - stop modulation')
        hold on
        xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
        colormap(jet)
        yline(length(find(speedModData{m}.nvm_stop_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        yline(length(find(speedModData{m}.nvm_stop_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
        
        set(gca, 'ydir', 'normal')
        
        speedModData{m}.pvm_stop_prePostDiff =  stop_prePostDiff(posCorr);
        speedModData{m}.nvm_stop_prePostDiff =  stop_prePostDiff(negCorr);
        
        subplot(3,4,4)
        hobj = bar(1,nanmean(speedModData{m}.pvm_stop_prePostDiff), 'FaceColor', [72,107, 142]./257);
        hold on
        errorbar(1, nanmean(speedModData{m}.pvm_stop_prePostDiff), nanstd(speedModData{m}.pvm_stop_prePostDiff)/length(speedModData{m}.pvm_stop_prePostDiff), 'k', 'linestyle', 'none', 'LineWidt', 2);
        hobj = bar(2, nanmean(speedModData{m}.nvm_stop_prePostDiff), 'FaceColor', [212, 147, 145]./257);
        errorbar(2, nanmean(speedModData{m}.nvm_stop_prePostDiff), nanstd(speedModData{m}.nvm_stop_prePostDiff)/length(speedModData{m}.nvm_stop_prePostDiff), 'k', 'linestyle', 'none', 'LineWidt', 2);
        xticks([1 2])
        xticklabels({'PVM', 'NVM'})
        set(gca, 'FontName', 'Gotham', 'FontSize', 18)
        ylabel('DFF_{post}-DFF_{pre}')
        title('running stop')
        ax = gca;
        ax.YAxis.Exponent = 0;
        
        colormap(jet)
        
    
     end
     
     if ~isempty(modulation)
         pvm_modulation  = modulation(posCorr); % correlation bootstrapping result
         nvm_modulation  = modulation(negCorr); % correlation bootstrapping result
         
         speedModData{m}.pvm_mod_run = modulation(posCorr);
         speedModData{m}.nvm_mod_run = modulation(negCorr);
         
         speedModData{m}.pvm_run_mod_val  = runningModulation(posCorr);
         speedModData{m}.nvm_run_mod_val  = runningModulation(negCorr);

         % sort them after modulation strength
         
         subplot(3,4,7)
         area(0:length(find(speedModData{m}.pvm_mod_run == -1)),[sort(speedModData{m}.pvm_run_mod_val(speedModData{m}.pvm_mod_run == -1)),0], 'FaceColor', [184, 90,51]./255, 'EdgeColor', [184, 90,51]./255)
         hold on
         if ~isempty(speedModData{m}.pvm_run_mod_val(speedModData{m}.pvm_mod_run == 0))
             area(length(find(speedModData{m}.pvm_mod_run == -1))+1:length(find(speedModData{m}.pvm_mod_run<= 0)), sort(speedModData{m}.pvm_run_mod_val(speedModData{m}.pvm_mod_run == 0)), 'FaceColor', [227,228,226]./255, 'EdgeColor', [227,228,226]./255)
         end
         if ~isempty(speedModData{m}.pvm_run_mod_val(speedModData{m}.pvm_mod_run == 1))
             area(length(find(speedModData{m}.pvm_mod_run <= 0)):length(speedModData{m}.pvm_mod_run)-1, sort(speedModData{m}.pvm_run_mod_val(speedModData{m}.pvm_mod_run == 1)), 'FaceColor', [201, 219, 170]./255, 'EdgeColor',[201, 219, 170]./255)
         end
         set(gca, 'FontName', 'Gotham', 'FontSize', 18)
         ylabel('DFF_{run}-DFF_{rest}')
         xlabel('# cell')
         title('mod. magnitude (PVM)')
         
         subplot(3,4,8)
         area(0:length(find(speedModData{m}.nvm_mod_run == -1)),[sort(speedModData{m}.nvm_run_mod_val(speedModData{m}.nvm_mod_run == -1)),0], 'FaceColor', [184, 90,51]./255, 'EdgeColor', [184, 90,51]./255)
         hold on
         area(length(find(speedModData{m}.nvm_mod_run == -1)):length(find(speedModData{m}.nvm_mod_run<= 0)),[ 0,sort(speedModData{m}.nvm_run_mod_val(speedModData{m}.nvm_mod_run == 0))], 'FaceColor', [227,228,226]./255, 'EdgeColor', [227,228,226]./255)
         if ~isempty(find(speedModData{m}.nvm_run_mod_val(speedModData{m}.nvm_mod_run == 1)))
             area(length(find(speedModData{m}.nvm_mod_run <= 0)):length(speedModData{m}.nvm_mod_run)-1, sort(speedModData{m}.nvm_run_mod_val(speedModData{m}.nvm_mod_run == 1)), 'FaceColor', [201, 219, 170]./255, 'EdgeColor',[201, 219, 170]./255)
         end
         set(gca, 'FontName', 'Gotham', 'FontSize', 18)
         ylabel('DFF_{run}-DFF_{rest}')
         xlabel('# cell')
         title('mod. magnitude (NVM)')
     end
    

    %% Reward effect calculation
    signalRew = mData(m).sData.imdata.roiSignals(2).dff;
    rewardInterval  = [];
    for k = 1:length(mData(m).sData.trials.trialStart)
        rewardInterval{k} = mData(m).sData.trials.trialStart(k)-62:mData(m).sData.trials.trialStart(k)+62;
    end
    
    rewSingleSignal = [];
    rewardSignal = [];
    reward_modulation_val = [];
    reward_modulation = [];
    for j = 1:size(signalRew,1)
        smoothSignal       = signalRew(j,:);%smoothdata(signal(j,:),2,'gaussian','SmoothingFactor', 0.85);
        for f = 1:length(mData(m).sData.trials.trialStart)
            rewSingleSignal(f,:) = smoothSignal(rewardInterval{f});
        end
        
        speedModData{m}.rewardSignal(j,:) = normalize(nanmean(rewSingleSignal,1),'range',[0,1]);
   
        reward_modulation_val(j) = nanmean(speedModData{m}.rewardSignal(j,63:93))- nanmean(speedModData{m}.rewardSignal(j,31:62));
        %% Bootstrapping
        allReward = reshape(cell2mat(rewardInterval)', 125,length(rewardInterval));
        pre_reward = allReward(31:62,:); all_pre_reward = pre_reward(:);
        post_reward = allReward(62:93,:); all_post_reward = post_reward(:);

        numBootstrap = 100;
        [~,bootsam_pre_reward]   = bootstrp(numBootstrap,@mean,all_pre_reward); %all_resting_indices_notReward);
        [~,bootsam_post_reward]  = bootstrp(numBootstrap,@mean,all_post_reward);

        %         signal = normalize(smoothSignal,'range',[0,1]);
        for i = 1:numBootstrap
            boot_pre_reward = all_pre_reward(bootsam_pre_reward(:,i));%all_resting_indices_notReward(bootsam_rest(:,i));
            boot_post_reward = all_post_reward(bootsam_post_reward(:,i));
            boot_reward_modulation(i) = nanmean(smoothSignal(boot_post_reward))- nanmean(smoothSignal(boot_pre_reward));

        end

        pd = fitdist(boot_reward_modulation','Normal');
        ci = paramci(pd,'Alpha',.01);    
        
        % A cell was considered to have a significant positive running modulation if 
        % CI is greater than 0, negative modulation if CI is less than 0, and no modulation if the CI contains 0.
         
        if ci(1,1) < 0 && ci(2,1) < 0
            reward_modulation(j) = -1;
        elseif ci(1,1)<= 0 && ci(2,1) >= 0
            reward_modulation(j) = 0;
        elseif ci(1,1) > 0 && ci(2,1) > 0
           reward_modulation(j) = 1;
        end

    end
    
   speedModData{m}.pvm_reward_modulation = reward_modulation(posCorr);speedModData{m}.pvm_reward_modulation_val = reward_modulation_val(posCorr);
   speedModData{m}.nvm_reward_modulation = reward_modulation(negCorr);  speedModData{m}.nvm_reward_modulation_val = reward_modulation_val(negCorr);
    
    speedModData{m}.pvm_rewSignal = speedModData{m}.rewardSignal(posCorr,:);    
    speedModData{m}.nvm_rewSignal = speedModData{m}.rewardSignal(negCorr,:);
    
    speedModData{m}.pvm_rewardCombined = [speedModData{m}.pvm_rewSignal(speedModData{m}.pvm_reward_modulation == 1,:);...
        speedModData{m}.pvm_rewSignal(speedModData{m}.pvm_reward_modulation == 0,:);...
        speedModData{m}.pvm_rewSignal(speedModData{m}.pvm_reward_modulation == -1,:)];
    
    imAlpha  = ones(size(speedModData{m}.pvm_rewardCombined));
    imAlpha(isnan(speedModData{m}.pvm_rewardCombined)) = 0;
    
    subplot(3,4,11)
    imagesc(speedModData{m}.pvm_rewardCombined, 'AlphaData', imAlpha)
    xticks([1 31 62 93 124])
    xticklabels([-2 -1 0 1 2])
    set(gca, 'FontName', 'Gotham', 'FontSize', 18)
    xlabel('time(s)')
    ylabel('cell')
    title('PVM - reward modulation')
    hold on
    xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
    colormap(jet)
    yline(length(find(speedModData{m}.pvm_reward_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
    yline(length(find(speedModData{m}.pvm_reward_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')

    speedModData{m}.nvm_rewardCombined = [speedModData{m}.nvm_rewSignal(speedModData{m}.nvm_reward_modulation == 1,:);...
        speedModData{m}.nvm_rewSignal(speedModData{m}.nvm_reward_modulation == 0,:);...
        speedModData{m}.nvm_rewSignal(speedModData{m}.nvm_reward_modulation == -1,:)];
    imAlpha  = ones(size(speedModData{m}.nvm_rewardCombined));
    imAlpha(isnan(speedModData{m}.nvm_rewardCombined)) = 0;
    
    subplot(3,4,12)
    imagesc(speedModData{m}.nvm_rewardCombined, 'AlphaData', imAlpha)
    xticks([1 31 62 93 124])
    xticklabels([-2 -1 0 1 2])
    set(gca, 'FontName', 'Gotham', 'FontSize', 18)
    xlabel('time(s)')
    ylabel('cell')
    title('NVM - reward modulation')
    hold on
    xline(62, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','--')
    colormap(jet)
    yline(length(find(speedModData{m}.nvm_reward_modulation == 1))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
    yline(length(find(speedModData{m}.nvm_reward_modulation >= 0))+0.5, 'LineWidth', 4, 'Color',[1 1 1], 'LineStyle','-')
    
%     saveas(fig,fullfile(savePathSpeed,'start_stop_sum'), 'fig')
%     saveas(fig,fullfile(savePathSpeed,'start_stop_sum'), 'png')
%     close(fig)
%     
    
    modPerc(1,:,m) = [length(find(speedModData{m}.pvm_start_modulation == 1)); length(find(speedModData{m}.pvm_start_modulation == 0)); length(find(speedModData{m}.pvm_start_modulation == -1))]./length(speedModData{m}.pvm_start_modulation);
    modPerc(2,:,m) = [length(find(speedModData{m}.pvm_stop_modulation == -1));length(find(speedModData{m}.pvm_stop_modulation == 0)); length(find(speedModData{m}.pvm_stop_modulation == 1))]./length(speedModData{m}.pvm_stop_modulation);
    modPerc(3,:,m) = [length(find(speedModData{m}.nvm_start_modulation == 1)); length(find(speedModData{m}.nvm_start_modulation == 0)); length(find(speedModData{m}.nvm_start_modulation == -1))]./length(speedModData{m}.nvm_start_modulation);
    modPerc(4,:,m) = [length(find(speedModData{m}.nvm_stop_modulation == -1)); length(find(speedModData{m}.nvm_stop_modulation == 0));length(find(speedModData{m}.nvm_stop_modulation == 1))]./length(speedModData{m}.nvm_stop_modulation);
    modPerc(5,:,m) = [length(find(speedModData{m}.pvm_reward_modulation == 1)); length(find(speedModData{m}.pvm_reward_modulation == 0)); length(find(speedModData{m}.pvm_reward_modulation == -1))]./length(speedModData{m}.pvm_reward_modulation);
    modPerc(6,:,m) = [length(find(speedModData{m}.nvm_reward_modulation == 1)); length(find(speedModData{m}.nvm_reward_modulation == 0)); length(find(speedModData{m}.nvm_reward_modulation == -1))]./length(speedModData{m}.nvm_reward_modulation);

    
   barfig = figure();
   b = bar(modPerc(:,:,m),'stacked');
   b(1).FaceColor = [204 219 175]./255; b(1).EdgeColor = 'none';
   b(2).FaceColor = [227 228 226]./255; b(2).EdgeColor = 'none';
   b(3).FaceColor = [172 96 61]./255; b(3).EdgeColor = 'none';
   
   set(gca, 'yticklabels', [])
   set(gca, 'xticklabels', [])
   close all
end




gol_pvm_pos_start_modulation = []; gol_pvm_no_start_modulation =[]; gol_pvm_neg_start_modulation = [];
gol_pvm_pos_stop_modulation = []; gol_pvm_no_stop_modulation =[]; gol_pvm_neg_stop_modulation = [];
gol_pvm_pos_rew_modulation = []; gol_pvm_no_rew_modulation =[]; gol_pvm_neg_rew_modulation = [];

gol_nvm_pos_start_modulation = []; gol_nvm_no_start_modulation =[]; gol_nvm_neg_start_modulation = [];
gol_nvm_pos_stop_modulation = []; gol_nvm_no_stop_modulation =[]; gol_nvm_neg_stop_modulation = [];
gol_nvm_pos_rew_modulation = []; gol_nvm_no_rew_modulation =[]; gol_nvm_neg_rew_modulation = [];

gol_pvm_run_mod_val = [];gol_nvm_run_mod_val = [];
gol_pvm_start_mod_val = [];gol_nvm_start_mod_val = [];
gol_pvm_stop_mod_val = []; gol_nvm_stop_mod_val = [];
gol_pvm_rew_mod_val = []; gol_nvm_stop_mod_val = [];


gol_pvm_reward_mod_val = [];
gol_nvm_reward_mod_val = [];


gol_pvm_pos_start_mod_val = []; gol_pvm_neg_start_mod_val = [];
gol_nvm_pos_start_mod_val = []; gol_nvm_neg_start_mod_val = [];
gol_pvm_pos_stop_mod_val = []; gol_pvm_neg_stop_mod_val = [];
gol_nvm_pos_stop_mod_val = []; gol_nvm_neg_stop_mod_val = [];

nvm_start  = [];
nvm_stop  = [];

pvm_start  = [];
pvm_stop  = [];
        
for i = 1:length(speedModData)
        gol_pvm_pos_start_modulation = [gol_pvm_pos_start_modulation; speedModData{i}.pvm_signal_in_start_interval(speedModData{i}.pvm_start_modulation == 1,:)];
        gol_pvm_no_start_modulation  = [gol_pvm_no_start_modulation; speedModData{i}.pvm_signal_in_start_interval(speedModData{i}.pvm_start_modulation == 0,:)];
        gol_pvm_neg_start_modulation = [gol_pvm_neg_start_modulation; speedModData{i}.pvm_signal_in_start_interval(speedModData{i}.pvm_start_modulation == -1,:)];

        gol_pvm_pos_stop_modulation = [gol_pvm_pos_stop_modulation; speedModData{i}.pvm_signal_in_stop_interval(speedModData{i}.pvm_stop_modulation == 1,:)];
        gol_pvm_no_stop_modulation = [gol_pvm_no_stop_modulation; speedModData{i}.pvm_signal_in_stop_interval(speedModData{i}.pvm_stop_modulation == 0,:)];
        gol_pvm_neg_stop_modulation = [gol_pvm_neg_stop_modulation; speedModData{i}.pvm_signal_in_stop_interval(speedModData{i}.pvm_stop_modulation == -1,:)];

        gol_pvm_pos_rew_modulation = [gol_pvm_pos_rew_modulation; speedModData{i}.pvm_rewSignal(speedModData{i}.pvm_reward_modulation == 1,:)];
        gol_pvm_no_rew_modulation = [gol_pvm_no_rew_modulation; speedModData{i}.pvm_rewSignal(speedModData{i}.pvm_reward_modulation == 0,:)];
        gol_pvm_neg_rew_modulation = [gol_pvm_neg_rew_modulation; speedModData{i}.pvm_rewSignal(speedModData{i}.pvm_reward_modulation == -1,:)];

        % nvm
        
        idxSession{i} = find(speedModData{i}.nvm_start_modulation == 1);
        gol_nvm_pos_start_modulation = [gol_nvm_pos_start_modulation; speedModData{i}.nvm_signal_in_start_interval(speedModData{i}.nvm_start_modulation == 1,:)];
        gol_nvm_no_start_modulation = [gol_nvm_no_start_modulation; speedModData{i}.nvm_signal_in_start_interval(speedModData{i}.nvm_start_modulation == 0,:)];
        gol_nvm_neg_start_modulation = [gol_nvm_neg_start_modulation; speedModData{i}.nvm_signal_in_start_interval(speedModData{i}.nvm_start_modulation == -1,:)];

        gol_nvm_pos_stop_modulation = [gol_nvm_pos_stop_modulation; speedModData{i}.nvm_signal_in_stop_interval(speedModData{i}.nvm_stop_modulation == 1,:)];
        gol_nvm_no_stop_modulation = [gol_nvm_no_stop_modulation; speedModData{i}.nvm_signal_in_stop_interval(speedModData{i}.nvm_stop_modulation == 0,:)];
        gol_nvm_neg_stop_modulation = [gol_nvm_neg_stop_modulation; speedModData{i}.nvm_signal_in_stop_interval(speedModData{i}.nvm_stop_modulation == -1,:)];

        gol_pvm_run_mod_val = [gol_pvm_run_mod_val speedModData{i}.pvm_run_mod_val];
        gol_nvm_run_mod_val = [gol_nvm_run_mod_val speedModData{i}.nvm_run_mod_val];

         gol_pvm_neg_start_mod_val = [gol_pvm_neg_start_mod_val speedModData{i}.pvm_start_prePostDiff(speedModData{i}.pvm_start_modulation == -1)];
         gol_nvm_pos_start_mod_val = [gol_nvm_pos_start_mod_val speedModData{i}.nvm_start_prePostDiff(speedModData{i}.nvm_start_modulation == 1)];%         %
        gol_pvm_pos_stop_mod_val = [gol_pvm_pos_stop_mod_val speedModData{i}.pvm_stop_prePostDiff(speedModData{i}.pvm_stop_modulation == 1)];
        gol_nvm_neg_stop_mod_val = [gol_nvm_neg_stop_mod_val speedModData{i}.nvm_stop_prePostDiff(speedModData{i}.nvm_stop_modulation == -1)];

        gol_pvm_start_mod_val(i)    = nanmean(speedModData{i}.pvm_start_prePostDiff);
        gol_nvm_start_mod_val(i)    = nanmean(speedModData{i}.nvm_start_prePostDiff);
        gol_pvm_stop_mod_val(i)     = nanmean(speedModData{i}.pvm_stop_prePostDiff);
        gol_nvm_stop_mod_val(i)     = nanmean(speedModData{i}.nvm_stop_prePostDiff);
        
        gol_pvm_reward_mod_val(i) = nanmean(speedModData{i}.pvm_reward_modulation_val);
        gol_nvm_reward_mod_val(i) = nanmean(speedModData{i}.nvm_reward_modulation_val);

        if ~isempty(speedModData{i}.speed_in_start_interval)
            meanspeed_in_start_interval(i,:) = nanmean(speedModData{i}.speed_in_start_interval,1);
            meanspeed_in_stop_interval(i,:)  = nanmean(speedModData{i}.speed_in_stop_interval,1);
        end

end

% 
% 
fig = figure('Position', [100 100 450 300]); 
xticks([0 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Speed (cm/s)')
xlim([0 124])
hold on

ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

F = 1:125;
amean = nanmean(meanspeed_in_start_interval,1);
astd = nanstd(meanspeed_in_start_interval,1)/sqrt(size(meanspeed_in_start_interval,1));
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0 0 0],'linestyle','none','FaceAlpha', 0.3);
hold on
plot(nanmean(meanspeed_in_start_interval,1),'k', 'LineWidth', 1.5)

ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

saveas(fig,fullfile(save_dir,'start_speed'), 'fig')
saveas(fig,fullfile(save_dir,'start_speed'), 'png')
saveas(fig,fullfile(save_dir,'start_speed'), 'pdf')


fig =figure('Position', [100 100 450 300]); 
xticks([0 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Speed (cm/s)')
xlim([0 124])
hold on
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

F = 1:125;
amean = nanmean(meanspeed_in_stop_interval,1);
astd = nanstd(meanspeed_in_stop_interval,1)/sqrt(size(meanspeed_in_stop_interval,1));;
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0 0 0],'linestyle','none','FaceAlpha', 0.3);
hold on

plot(nanmean(meanspeed_in_stop_interval,1),'k', 'LineWidth', 1.5)

ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

saveas(fig,fullfile(save_dir,'stop_speed'), 'fig')
saveas(fig,fullfile(save_dir,'stop_speed'), 'png')
saveas(fig,fullfile(save_dir,'stop_speed'), 'pdf')

%% start figure: 
fig = figure('Position', [100 100 450 300]); 

[~, sortIdxPos] = sort(max(gol_pvm_pos_start_modulation'));
[~, sortIdxNeg] = sort(max(gol_pvm_neg_start_modulation'));
[~, sortIdxNo] = sort(max(gol_pvm_no_start_modulation'));

matrPos = [gol_pvm_pos_start_modulation(sortIdxPos,:);gol_pvm_no_start_modulation(sortIdxNo,:) ; gol_pvm_neg_start_modulation(sortIdxNeg,:) ];


imagesc(matrPos )
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Cell')
title('PM - start modulation')
colormap(jet)
hold on
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
colorbar

saveas(fig,fullfile(save_dir,'PM_start'), 'fig')
saveas(fig,fullfile(save_dir,'PM_start'), 'png')
saveas(fig,fullfile(save_dir,'PM_start'), 'pdf')


fig =figure('Position', [100 100 450 300]); 

[~, sortIdxPos] = sort(max(gol_nvm_pos_start_modulation'));
[~, sortIdxNeg] = sort(max(gol_nvm_neg_start_modulation'));
[~, sortIdxNo] = sort(max(gol_nvm_no_start_modulation'));

matrNeg  = [gol_nvm_pos_start_modulation(sortIdxPos,:); gol_nvm_no_start_modulation(sortIdxNo,:); gol_nvm_neg_start_modulation(sortIdxNeg,:) ];

imagesc(matrNeg)
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Cell')
title('NM - start modulation')
colormap(jet)
hold on
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

saveas(fig,fullfile(save_dir,'NM_start'), 'fig')
saveas(fig,fullfile(save_dir,'NM_start'), 'png')
saveas(fig,fullfile(save_dir,'NM_start'), 'pdf')
colorbar

fig =figure('Position', [100 100 450 300]); 

[~, sortIdxPos] = sort(max(gol_pvm_pos_stop_modulation'));
[~, sortIdxNeg] = sort(max(gol_pvm_neg_stop_modulation'));
[~, sortIdxNo] = sort(max(gol_pvm_no_stop_modulation'));

matrPos = [gol_pvm_pos_stop_modulation(sortIdxPos,:);gol_pvm_no_stop_modulation(sortIdxNo,:) ; gol_pvm_neg_stop_modulation(sortIdxNeg,:) ];


imagesc(matrPos )
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Cell')
title('PM - stop modulation')
colormap(jet)
hold on
colorbar
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

saveas(fig,fullfile(save_dir,'PM_stop'), 'fig')
saveas(fig,fullfile(save_dir,'PM_stop'), 'png')
saveas(fig,fullfile(save_dir,'PM_stop'), 'pdf')

fig =figure('Position', [100 100 450 300]); 

[~, sortIdxPos] = sort(max(gol_nvm_pos_stop_modulation'));
[~, sortIdxNeg] = sort(max(gol_nvm_neg_stop_modulation'));
[~, sortIdxNo] = sort(max(gol_nvm_no_stop_modulation'));

matrNeg  = [gol_nvm_pos_stop_modulation(sortIdxPos,:); gol_nvm_no_stop_modulation(sortIdxNo,:); gol_nvm_neg_stop_modulation(sortIdxNeg,:) ];

imagesc(matrNeg)
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time (s)')
ylabel('Cell')
title('NM - stop modulation')
colormap(jet)
hold on
colorbar
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

saveas(fig,fullfile(save_dir,'NM_stop'), 'fig')
saveas(fig,fullfile(save_dir,'NM_stop'), 'png')
saveas(fig,fullfile(save_dir,'NM_stop'), 'pdf')



figure('Position', [100 100 200 300]); 

subplot(2,1,1)
hobj = bar(1,nanmean(gol_pvm_start_mod_val), 'FaceColor', [72,107, 142]./257);
hold on
errorbar(1, nanmean(gol_pvm_start_mod_val), nanstd(gol_pvm_start_mod_val)/sqrt(length(speedModData)), 'k', 'linestyle', 'none', 'LineWidt', 2);
hobj = bar(2, nanmean(gol_nvm_start_mod_val), 'FaceColor',[68 195 214]./257 );
errorbar(2, nanmean(gol_nvm_start_mod_val), nanstd(gol_nvm_start_mod_val)/sqrt(length(speedModData)), 'k', 'linestyle', 'none', 'LineWidt', 2);
xticks([1 2])
xticklabels({'PVM', 'NVM'})
set(gca, 'FontName', 'Arial', 'FontSize', 14)
ylabel('\Delta F/F_{post}-\Delta F/F_{pre}')
title('running start')
ax = gca;
ax.YAxis.Exponent = 0;
hold on
xlim([0 3])    
box off

ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];

subplot(2,1,2)
hobj = bar(1,nanmean(gol_pvm_stop_mod_val), 'FaceColor', [72,107, 142]./257);
hold on
errorbar(1, nanmean(gol_pvm_stop_mod_val), nanstd(gol_pvm_stop_mod_val)/sqrt(length(speedModData)), 'k', 'linestyle', 'none', 'LineWidt', 2);
hobj = bar(2, nanmean(gol_nvm_stop_mod_val), 'FaceColor', [68 195 214]./257 );
errorbar(2, nanmean(gol_nvm_stop_mod_val), nanstd(gol_nvm_stop_mod_val)/sqrt(length(speedModData)), 'k', 'linestyle', 'none', 'LineWidt', 2);
xticks([1 2])
xticklabels({'PVM', 'NVM'})
set(gca, 'FontName', 'Arial', 'FontSize', 14)
ylabel('\Delta F/F_{post}-\Delta F/F_{pre}')
title('running stop')
ax = gca;
ax.YAxis.Exponent = 0;
hold on
    
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
xlim([0 3])    
box off

