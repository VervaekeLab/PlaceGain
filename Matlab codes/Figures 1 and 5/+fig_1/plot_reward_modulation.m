rng(1)
save_dir = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/analysis/VIP-GCAMP/lag_1s';

% load(fullfile(save_dir,'neg_chance'))
% load(fullfile(save_dir,'pos_chance'))

speedLim = 10; 
slowIntervalLength = 31;
pvm_rewSignal = [];
nvm_rewSignal = [];
pvm_pos_rewSignal =[];
pvm_neg_rewSignal =[];
nvm_pos_rewSignal =[];
nvm_neg_rewSignal =[];
speed_all = [];
water_all = [];
licking_all = [];
speedModData = cell(length(sessions));
for m = 1:length(sessions)
    smoothedPlaceCellSignal = [];
    signal      = mData(m).sData.imdata.roiSignals(2).dff;
    speed       = mData(m).sData.behavior.runSpeedDs;
    waterReward = mData(m).sData.behavior.waterRewardDs;
    licking     = mData(m).sData.behavior.lickDs;
    

    signal(deletedROIs{m},:)   = [];
    minLength = min([size(signal,2),length(speed),length(waterReward)]);

    signal = signal(:,1:minLength);
    speed = speed(1:minLength);
    waterReward = waterReward(1:minLength);
    position = mData(m).sData.behavior.wheelPosDs(2:minLength+1); 
    
   
    %% NOW: Bootstrapping to check significance for stop and start modulation
    
      
    % load correlation values
    load(strcat(save_dir,'/',sessions{m}(1:5),'/',sessions{m},'/speedCorrSum.mat'))
    
%     pos_chance = 0.1875;
%     neg_chance =-0.1865;
    posCorr = find(speedCorr{m}.val > pos_chance);
    negCorr = find(speedCorr{m}.val < neg_chance);
    

    %% Reward effect calculation
    signalRew = mData(m).sData.imdata.roiSignals(2).dff;
    rewardInterval  = [];
    for k = 1:length(mData(m).sData.trials.trialStart)
        rewardInterval{k} = mData(m).sData.trials.trialStart(k)-62:mData(m).sData.trials.trialStart(k)+62;
        speed_during_reward(k,:) = speed(rewardInterval{k});
        water_during_reward(k,:) = waterReward(rewardInterval{k});
        licking_during_reward(k,:) = licking(rewardInterval{k});
    end
    
    rewSingleSignal = [];
    rewardSignal = [];
    reward_modulation_val = [];
    reward_modulation = [];
    reward_modulation_bootstrapped=[];
    for j = 1:size(signalRew,1)
        smoothSignal       = signalRew(j,:);%smoothdata(signal(j,:),2,'gaussian','SmoothingFactor', 0.85);
        for f = 1:length(mData(m).sData.trials.trialStart)
            rewSingleSignal(f,:) = smoothSignal(rewardInterval{f});
        end
        
        speedModData{m}.rewardSignal(j,:) = normalize(nanmean(rewSingleSignal,1),'range',[0,1]);
   
        reward_modulation_val(j) = nanmean(speedModData{m}.rewardSignal(j,32:93))- nanmean([speedModData{m}.rewardSignal(j,1:31),speedModData{m}.rewardSignal(j,94:end)]);
        %% Bootstrapping
        allReward = reshape(cell2mat(rewardInterval)', 125,length(rewardInterval));
        during_reward = allReward(32:93,:); all_pre_reward = pre_reward(:);
        outside_reward =  [allReward(1:31,:);allReward(94:end,:)]; all_post_reward = outside_reward(:);
        [speedModData{m}.rewardSignal(j,1:31),speedModData{m}.rewardSignal(j,93:end)];
        all_post_reward = post_reward(:);
        allReward(32:62,:); all_post_reward = post_reward(:);

        numBootstrap = 3000;
        [~,bootsam_pre_reward]   = bootstrp(numBootstrap,@mean,during_reward); %all_resting_indices_notReward);
        [~,bootsam_post_reward]  = bootstrp(numBootstrap,@mean,outside_reward);

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
            reward_modulation_bootstrapped(j) = -1;
        elseif ci(1,1)<= 0 && ci(2,1) >= 0
            reward_modulation_bootstrapped(j) = 0;
        elseif ci(1,1) > 0 && ci(2,1) > 0
           reward_modulation_bootstrapped(j) = 1;
        end

    end
    
    reward_modulation(reward_modulation_val>0) = 1;
    reward_modulation(reward_modulation_val<0) = -1;
    
    reward_final = reward_modulation'+reward_modulation_bootstrapped';

    pos_reward_signal = speedModData{m}.rewardSignal(posCorr,:);
    neg_reward_signal = speedModData{m}.rewardSignal(negCorr,:);

    pvm_pos_rewSignal = [pvm_pos_rewSignal;pos_reward_signal(reward_final(posCorr)==2,:)];
    nvm_pos_rewSignal = [nvm_pos_rewSignal;neg_reward_signal(reward_final(negCorr)==2,:)];
    
    pvm_neg_rewSignal = [pvm_neg_rewSignal;pos_reward_signal(reward_final(posCorr)==-2,:)];
    nvm_neg_rewSignal = [nvm_neg_rewSignal;neg_reward_signal(reward_final(negCorr)==-2,:)];
    
    speed_all = [speed_all;nanmean(speed_during_reward,1)];
    water_all = [water_all;nanmean(water_during_reward,1)];
    licking_all = [licking_all;nanmean(licking_during_reward,1)];

end
% 
% imAlpha  = ones(size(cell2mat(speedModData{m}.pvm_rewSignal)));
% imAlpha(isnan(cell2mat(speedModData{m}.pvm_rewSignal))) = 0;
[~, sortIdxPVM_pos] = sort(max(pvm_pos_rewSignal'));
[~, sortIdxNVM_pos] = sort(max(nvm_pos_rewSignal'));
[~, sortIdxPVM_neg] = sort(max(pvm_neg_rewSignal'));
[~, sortIdxNVM_neg] = sort(max(nvm_neg_rewSignal'));

figure()
subplot(2,2,1)
imagesc(pvm_pos_rewSignal(sortIdxPVM_pos,:))
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
set(gca, 'FontName', 'Arial', 'FontSize', 16)
xlabel('Time(s)')
ylabel('Cell #')
title('PVM - positive modulation')
hold on

subplot(2,2,2)
imagesc(pvm_neg_rewSignal(sortIdxPVM_neg,:))
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
set(gca, 'FontName', 'Arial', 'FontSize', 16)
xlabel('Time(s)')
ylabel('Cell #')
title('PVM - negative modulation')
hold on

subplot(2,2,3)
imagesc(nvm_pos_rewSignal(sortIdxNVM_pos,:))
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
set(gca, 'FontName', 'Arial', 'FontSize', 16)
xlabel('Time(s)')
ylabel('Cell #')
title('NVM - positive modulation')
hold on

subplot(2,2,4)
imagesc(nvm_neg_rewSignal(sortIdxNVM_neg,:))
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
set(gca, 'FontName', 'Arial', 'FontSize', 16)
xlabel('Time(s)')
ylabel('Cell #')
title('NVM - negative modulation')
hold on


figure()
subplot(3,1,1)
F = 1:125;
amean = nanmean(speed_all,1);
astd = nanstd(speed_all,1);
plot(amean,'k','LineWidth',1.5)
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0 0 0],'linestyle','none','FaceAlpha', 0.3);
xlim([0 124])
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
xlabel('Time(s)')
box off
ylabel('Speed (cm/s)')
set(gca, 'FontName', 'Arial', 'FontSize', 16)
subplot(3,1,2)
amean = nanmean(water_all,1);
astd = nanstd(water_all,1);
plot(amean,'k','LineWidth',1.5)
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0 0 0],'linestyle','none','FaceAlpha', 0.3);
box off
% plot(nanmean(water_all),'k','LineWidth',1.5)
xlim([0 124])
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
yticks([0 1])
ylabel('Water delivery(on/off)')

xlabel('Time(s)')
set(gca, 'FontName', 'Arial', 'FontSize', 16)
subplot(3,1,3)
amean = nanmean(licking_all,1);
astd = nanstd(licking_all,1);
plot(amean,'k','LineWidth',1.5)
patch([F(~isnan(amean)) fliplr(F(~isnan(amean)))],[amean(~isnan(amean))+astd(~isnan(amean)) fliplr(amean(~isnan(amean))-astd(~isnan(amean)))],[0 0 0],'linestyle','none','FaceAlpha', 0.3);
xlim([0 124])
xticks([1 31 62 93 124])
xticklabels([-2 -1 0 1 2])
box off
yticks([0:0.2:0.8])
yticks(100*[0:0.2:0.8])
ylabel('Lick probability (%)')
xlabel('Time(s)')
set(gca, 'FontName', 'Arial', 'FontSize', 16)


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

