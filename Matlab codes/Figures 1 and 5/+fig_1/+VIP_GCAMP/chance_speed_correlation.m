
% set fixed seed
rng(1);

n_shuffle     = 1000;
lag_time_bins = 31;

for m = 1:length(sessions)
    smoothedPlaceCellSignal = [];
   
    signal = mData(m).sData.imdata.roiSignals(2).dff;
    speed = smoothdata(normalize(mData(m).sData.behavior.runSpeedDs,'range',[0 1]),'gaussian',10);
    
    % sometimes the speed is 1 or two indices longer
    min_length = min([length(speed),size(signal,2)]);
    speed = speed(1:min_length);
    signal = signal(:,1:min_length);

    signal(:,isnan(speed)) = [];
    speed(isnan(speed)) = [];

    signal(deletedROIs{m},:)   = [];
 
    speedCorr_chance{m}.val = [];
    speedCorr_chance{m}.c = [];
    speedCorr_chance{m}.idx  = [];

    for j = 1: size(signal,1)
        f = 0;
        for k = 1:n_shuffle
            f = f+1;
            smoothSignal       = circshift(smoothdata(signal(j,:),'gaussian',10),randi(size(signal,2),1));%smoothdata(signal(j,:),'gaussian','SmoothingFactor', 0.85),'range',[0,1]);
            [speedCorr_chance{m}.c(:,j),speedCorr_chance{m}.lags] = xcorr(speed-mean(speed),smoothSignal'-mean(smoothSignal), lag_time_bins,'coeff');
            [~ ,speedCorr_chance{m}.idx(f)] = max(abs(speedCorr_chance{m}.c(:,j)));
            speedCorr_chance{m}.val(f) = speedCorr_chance{m}.c(speedCorr_chance{m}.idx(j),j);
        end
    end

    speedCorr_chance{m}.posCorr = find(speedCorr_chance{m}.val(:)> 0);
    speedCorr_chance{m}.negCorr = find(speedCorr_chance{m}.val(:) < 0);
   
end


corr_Pos = []; corr_Neg = []; 
for i = 1:length(sessions)
    corr_Pos = [corr_Pos, speedCorr_chance{i}.val(speedCorr_chance{i}.posCorr)];
    corr_Neg = [corr_Neg, speedCorr_chance{i}.val(speedCorr_chance{i}.negCorr)];
end
    


fig = figure();
hx = histogram(corr_Pos,'Facecolor',[72,107, 142]./257, 'EdgeColor', 'none');
hold on
bin_Width = hx.BinWidth;
histogram(corr_Neg, 'Facecolor',[212, 147, 145]./257, 'EdgeColor', 'none', 'BinWidth', bin_Width);
set(gca, 'FontName', 'Gotham', 'FontSize',15)
xlabel('Peak velocity cross correlation')
ylabel('# cells')

pd = fitdist([corr_Pos,corr_Neg]','Normal');
civals = paramci(pd);
xline(nanmean([corr_Pos,corr_Neg])-2*nanstd([corr_Pos,corr_Neg]))
xline(nanmean([corr_Pos,corr_Neg])+2*nanstd([corr_Pos,corr_Neg]))
neg_chance = nanmean([corr_Pos,corr_Neg])-2*nanstd([corr_Pos,corr_Neg]);
pos_chance = nanmean([corr_Pos,corr_Neg])+2*nanstd([corr_Pos,corr_Neg]);

if ~exist(fullfile(save_dir,'speedCorr_chance'))
    mkdir(fullfile(save_dir,'speedCorr_chance'));
end
saveas(fig,fullfile(save_dir,'speedCorr_chance'), 'fig')
saveas(fig,fullfile(save_dir,'speedCorr_chance'), 'png')
saveas(fig,fullfile(save_dir,'speedCorr_chance'), 'pdf')

save(fullfile(save_dir,'speedCorr_chance','neg_chance.mat'),'neg_chance')
save(fullfile(save_dir,'speedCorr_chance','pos_chance.mat'),'pos_chance')

% 