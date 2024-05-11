

file =  { 'm8137-20240117-01','m8138-20240117-01','m8139-20240117-02',...
                                'm8140-20240117-01'};

dirct ='/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/decoding/VIP-Chr-HPC/';

% we used this to discard sessions in which the control decoding value was
% not better than the 5th quantile of the shuffled distribution

rng(1)
numShuff = 1000;

for i =1:length(file)
    
    
    directory = fullfile(dirct,file{i});
    
    load(strcat(directory,'/d_data_off.mat'));
    
    rmseError = [];

    for m = 1:5
        realTime  = d_data.iter{d_data.medianTestFold}.realPos_test;
        predTime   = d_data.iter{d_data.medianTestFold}.predPos_test;
        
        for f = 1:numShuff
            circPos        = circshift(realTime  ,randi(length(realTime  ),1));
            cirPred        = circshift(predTime  ,randi(length(predTime  ),1));
            rmseError(m,f)= calculateTimeError(circPos,predTime);
        end
    end
    
    nanmean(d_data.rmse_test)
    mean(rmseError(:))
    quantile(rmseError(:),0.05)
    
end
 
    




function mean_predError = calculateTimeError(actualClass,predClass)
% unique(actualClass)

maxval= max(actualClass);
for i = 1:length(predClass)
    predError(i) = abs(actualClass(i)-predClass(i));
    if predError(i) > maxval/2 && actualClass(i) > predClass(i)
        predError(i) = abs(maxval-actualClass(i)+predClass(i));
    elseif predError(i) > maxval/2 && actualClass(i) < predClass(i)
        predError(i) = abs(maxval-predClass(i)+actualClass(i));
    end
end

mean_predError = nanmean(predError);
sem_predError  = nanstd(predError)/sqrt(length(predError));
end
