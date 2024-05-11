clearvars -except p files file load_dirct save_dirct  sessionInfo
close all

 % File paths
dataDirectory = fullfile(load_dirct, file(1:5), file, '/');
saveDirct = fullfile(save_dirct, file(1:5), file);

% Create directory if it doesn't exist
if ~exist(saveDirct, 'dir')
    mkdir(saveDirct);
end

% Load data
load(fullfile(dataDirectory, 'sData.mat'));

% Parameters
params.deletePosThreshold = 20;
params.binSize = 10.4667;
params.timebinSize = 1;
params.smooth = 0;
params.intensity_idx = 1; % in some sessions,there werre more than one intensity used, replace intensity_idx
params.numFolds = 5;  % set number of cross validation folds
params.binGaps  = params.binSize/2:params.binSize:157.5; % create evenly sized bins, that consider circularity

% set no of cells you want to use for decoding:
noCells = size(sData.imdata.roiSignals(2).dff,1);
params.cellIdx = 1:noCells; % used cells are randomly chosen from full set
signal = normalize(double(sData.imdata.roiSignals(2).dff(params.cellIdx,:)),1 ,'range', [0 1]);
position    = sData.behavior.wheelPosDs;

% Extract signalOpto
signalOpto = sData.daqdata.optoSignal(sData.daqdata.frameIndex);

%% Provide trials
nTrials     = numel(sData.trials.trialStart)-1;
optoTrials  = sData.behavior.opto.OptoOnTrials;%(length(sData.behavior.opto.OptoOnTrials)-nTrials+1:end);
afterTrials = sData.behavior.opto.AfterOptoTrials;%(length(sData.behavior.opto.OptoOnTrials)-nTrials+1:end);
offTrials   = sData.behavior.opto.OptoOffTrials;%(length(sData.behavior.opto.OptoOnTrials)-nTrials+1:end);

% Find trial indices
offTrialsIdx    = find(offTrials==1);
afterTrialsIdx  = find(afterTrials == 1);
optoTrialsIdx   = find(optoTrials == 1);

% in some sessions,there werre more than one intensity used, replace intensity_idx
% optoTrialsIdx = find(any(sData.behavior.opto.OptoStimOnMatrix(end - nTrials + 1:end, :) == params.intensity_idx, 2));

% make sure that the optogenetic laser was defineylt not on during off
% trials
for i = 1:length(offTrialsIdx)
    isInTrial   = sData.trials.trialLength == offTrialsIdx(i);
    if nanmean(signalOpto(isInTrial)) > 1
        offTrials(offTrialsIdx(i)) = NaN;
    end
end
offTrialsIdx = find(offTrials==1);   


minTrials       = min([length(optoTrialsIdx) length(offTrialsIdx) length(afterTrialsIdx)]);
optoTrialsIdx   = sort(optoTrialsIdx(randperm(length(optoTrialsIdx),minTrials)));
offTrialsIdx    = sort(offTrialsIdx(randperm(length(offTrialsIdx),minTrials)));
afterTrialsIdx  = sort(afterTrialsIdx(randperm(length(afterTrialsIdx),minTrials)));


minLength = min([size(signal,2),length(sData.trials.trialLength)]);
sData.trials.trialLength =sData.trials.trialLength(1:minLength);
% find indices for each experimental state
offIdx = [];   for i = 1:length(offTrialsIdx);   offIdx =[offIdx find(sData.trials.trialLength == offTrialsIdx(i))]; end
afterIdx = []; for i = 1:length(afterTrialsIdx); afterIdx =[afterIdx find(sData.trials.trialLength == afterTrialsIdx(i))]; end
optoIdx = [];  for i = 1:length(optoTrialsIdx);  optoIdx =[optoIdx find(sData.trials.trialLength == optoTrialsIdx(i))]; end


d_data_opto = [];
d_data_after = [];
d_data_off = [];

%%
optoSignal  = signal(:,optoIdx);
afterSignal = signal(:,afterIdx);
offSignal   = signal(:,offIdx);

optoPos  = position(optoIdx);  
afterPos = position(afterIdx);
offPos   = position(offIdx);  

% Bin the position and signals
[afterPosBinned, afterSignalBinned,params] = binData(afterPos, afterSignal, params);
[optoPosBinned, optoSignalBinned,params]   = binData(optoPos, optoSignal, params);
[offPosBinned, offSignalBinned,params]     = binData(offPos, offSignal, params);

% divide the data up into 5*num_folds pieces
bins = 1:length(params.binGaps);

saveDecodingData(afterSignalBinned, afterPosBinned,params, saveDirct, '/d_data_after.mat');
saveDecodingData(offSignalBinned, offPosBinned, params,saveDirct, '/d_data_off.mat');
saveDecodingData(optoSignalBinned, optoPosBinned,params, saveDirct, '/d_data_opto.mat');


function data = findEmptyROIs(data)

deleteCell = [];
for n = 1:length(data.cellIdx)
   if length(find(isnan(data.deconvTrain(n,:)) == 1)) == length(data.deconvTrain(n,:)) ||length(find(isnan(data.deconvTest(n,:)) == 1)) == length(data.deconvTest(n,:)) 
        deleteCell = [deleteCell; n];
    end
%     if length(find(isnan(data.deconvTrain(n,:)) == 1)) > 10 ||length(find(isnan(data.deconvTest(n,:)) == 1)) > 10
%         deleteCell = [deleteCell; n];
%     end
end

data.deconvTrain(deleteCell,:) = []; data.deconvTest(deleteCell,:) = [];
data.cellIdx(deleteCell) = [];

end

function [data deleteIdxTrain deleteIdxTest] = findEmptyIdx(data)

deleteIdxTrain = [];
for n = 1:size(data.deconvTrain,2)
    if length(find(isnan(data.deconvTrain(:,n)) == 1)) > 50 || nanmean(data.deconvTrain(:,n))> 5*nanmean(data.deconvTrain(:))
        deleteIdxTrain = [deleteIdxTrain; n];
    
    end
end

data.deconvTrain(:,deleteIdxTrain) = []; 

deleteIdxTest = [];
for n = 1:size(data.deconvTest,2)
    if length(find(isnan(data.deconvTest(:,n)) == 1)) > 50 || nanmean(data.deconvTest(:,n))> 5*nanmean(data.deconvTest(:))
        deleteIdxTest = [deleteIdxTest; n];
    end
end

data.deconvTest(:,deleteIdxTest) = [];
% data.cellIdx(deleteCell) = [];

end

function lambda = poissonGLMPos(signal, pos, smooth)

% w = linspace(0,1,16);
% w = flipud(w(2:end)');

% f = zeros(size(signal,1),max(pos)+1); % 50x78
% lambda = zeros(size(signal,1),max(pos)); % 59x 77
f = zeros(size(signal,1),max(pos)+1);
% f= zeros(size(signal,1), 76);
X = zeros(length(pos),max(pos)+1);
lambda = zeros(size(signal,1),max(pos));

for i = 1:size(signal,1)
        X = zeros(length(pos),length(unique(pos)+1));
        for j = 1:length(pos); X(j,pos(j)+1) = 1; end;  X(:,1) = 1;
        
        
        % find optimising parameters
        init         = (X'*signal(i,:)')/sum(signal(i,:));init(isnan(init)) = 0;
        opts         = optimset('Gradobj', 'on','Hessian','off','Display','off');
        lossfun      = @(params) neglogli_poissGLM(params, X, signal(i,:)', smooth);
        f(i,:)       = fminunc(lossfun, init, opts);
        
        c                = f(i,1);
        param            = f(i,2:end);
        X                = eye(size(X,2)-1,size(X,2)-1);        
        lambda(i,:)      = exp(c +X*param');
       
end

end

function [neglogli, dL, H] = neglogli_poissGLM(params,X,Y, beta)
% Compute negative log-likelihood of data undr Poisson GLM model with
% exponential nonlinearity

vv = X*params; 
rr = exp(vv);%*dtbin;

% beta = 1;
[J,J_g,J_h] = penalty(params,beta);
% ------------  Compute log-likelihood -----------
term1 = -vv'*Y;
term2 = sum(rr);  
neglogli = term1+term2+ J ;

% % ---------  Compute Gradient - First Derivative --
if (nargout > 1)
    dL1 = -X'*Y; 
    dL0 = X'*rr; 
    dL = dL1+dL0+ J_g;    
end

% [J,J_g,J_h] = penalty(param,b);

%% compute f, the gradient, and the hessian 
% f = sum(rate-Y(:).*u) + J ;
% df = real(X' * (rate - Y(:)) + J_g);
% hessian = hessian_glm + blkdiag(J_h);
% 


% ----------  Compute Hessian - Second Derivative ----
if nargout > 2
    H = XX'*bsxfun(@times,X,rr)+blkdiag(J_h); 
end

end

function [J,J_g,J_h] = penalty(param,beta)
    %% params.smoothing functions called in the above script
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    J = beta*0.5*param'*DD1*param;
    J_g = beta*DD1*param;
    J_h = beta*DD1;
end

function likelihood= llh_signalGivPos(deconv,f)
%lambda1(n,:) = b(n)+g(n)*exp(k(n)*cos(x-mu(n)));
% calculate  P(spike(i)| time(i)) 

prob = zeros(size(f)); 
likelihood = zeros(size(deconv,2),size(f,2));

for t = 1: size(deconv,2) 
    % select signal for one time bin
    frTemp = repmat(deconv(:,t),1,size(f,2));
%     for j = 1:size(frTemp,2)
%         for i = 1:size(frTemp,1)
%             %prob. of the time points signals (separate per neuron) given position
%             prob(i,j) = exp(frTemp(i,j) * log(f(i,j)) - f(i,j) - gammaln(frTemp(i,j)+1)); 
%         end
%     end
    
    prob = exp(frTemp .* log(f) - f - gammaln(frTemp+1)); 

    % probability signal(each time point) given position, multiply over
    % neurons -> here we assume statistical independence
    likelihood(t,:) = prod(prob); %Prob. of all signals (all neurons) given pos
end

% normalisation

likelihood = likelihood./nansum(likelihood,2);

end

function classificationAccuracy = confusionMatrix(predClass, actualClass, directory, filename, classes)

% classes = unique(actualClass);
confMatrix = zeros(length(classes),length(classes));
corrPred = zeros(length(classes),1);

for i = 1: length(classes)
    idx                   = actualClass == classes(i);
    predClasstemp         = predClass(idx);
    corrPred(i)           = numel(find(predClasstemp == classes(i)));
    
    confMatrix(i,i)       = corrPred(i);
    
    falsePred = zeros(length(classes),1);
    for j = 1: length(classes)
        if j == i
        else
            falsePred(j) = length(find(predClasstemp == classes(j)));
            confMatrix(j,i) = falsePred(j);
        end
    end
end

classificationAccuracy = sum(corrPred) / length(predClass);

fig          = figure;
fig.Name     = 'Confusion Matrix';

% subplot(2,1,1)
imagesc(flipud(confMatrix));
xlabel('actual class');
ylabel('predicted class');
% set(gca,'yTick',[1 11 21 31], 'yTickLabel',{'30','20', '10', '0'});
set(gca,'FontSize', 18, 'FontName','Gotham');
colorbar()

% predPos = zeros(length(classes),1);
% error = zeros(length(classes),1);
% for i = 1: length(classes)
%     idx = find(actualClass == classes(i));
%     predClassIdx = predClass(idx);
%     if classes(i) < mean(classes)
%         idxChange1 = find(predClassIdx-classes(i)> mean(classes));
%         predClassIdx(idxChange1) = max(classes)-predClassIdx(idxChange1);
%     elseif classes(i) > mean(classes)
%         idxChange2 = find(predClassIdx-classes(i) < -mean(classes));
%         predClassIdx(idxChange2) = classes(i)+predClassIdx(idxChange2);
%     end
%     predPos(i) = nanmean(predClassIdx);
%     error(i)   = nanstd(predClassIdx);
% end

% subplot(2,1,2)
% errorbar(unique(classes), predPos, error);
% xlabel('actual position (cm)');
% ylabel('predicted position (cm)');
% set(gca,'FontName','Gotham','FontSize',18);
% set(gca, 'yTick', [0 10 20 30 40 50], 'yTickLabel', {'0','40','80', '120','160'},...
%     'xTick', [0 10 20 30 40 50], 'xTickLabel', {'0','40','80', '120','160'})
% hold on
% plot(unique(classes), unique(classes), '--')
% ylim([min(classes), max(classes)]);

saveas(fig,fullfile(directory,filename))
saveas(fig,fullfile(directory,filename),'png')
close(fig)

end

function rmse = plotDecodingPos(predPos, pos, keepClose, directory, filename, maxval)

fig = figure();
subplot(2,1,1);
plot(1:length(pos),pos)
hold on
plot(1:length(predPos),predPos, 'LineWidth', 1, 'Color', [198, 82, 77]./255)
legend('real position','predicted position', 'eastoutside');

xlim([1 length(predPos)])
xlabel('Pos(Samples)')
ylabel('Position in [cm]')
yticks([0 40 80 120 160 ])
yticklabels({'0','40','80','120','160'})
set(gca, 'FontName', 'Gotham');
set(gca, 'FontSize' ,18);

predError = zeros(length(pos),1);
for i = 1:length(pos)
    predError(i) = abs(pos(i)-predPos(i));
    if predError(i) > maxval/2 && pos(i) > predPos(i)
        predError(i) = abs(maxval-pos(i)+predPos(i));
    elseif predError(i) > maxval/2 && pos(i) < predPos(i)
        predError(i) = abs(maxval-predPos(i)+pos(i));
    end
end

rmse = mean(predError);

subplot(2,1,2);
plot(1:length(pos),predError, 'LineWidth',1, 'Color', [49, 54, 83]./255)

xlim([1 length(predPos)])
xlabel('Pos(Samples)')
ylabel('Error in [cm]')
yticks([0 40 80])
yticklabels({'0','40','80'})
set(gca, 'FontName', 'Gotham');
set(gca, 'FontSize' ,18);

if strcmp(keepClose,'close') == 1
    saveas(fig,fullfile(directory,filename))
    close
elseif strcmp(keepClose, 'closeNotSave') == 1
    close
end

end

function d_data = main_decoding(deconvSignal,posBinned,params)

    sections = params.numFolds*5;
    edges = round(linspace(1,size(deconvSignal,2)+1,sections+1));

    rmse_train = zeros(1,params.numFolds);
    rmse_test = zeros(1,params.numFolds);
    d_data = [];
    
    for i = 1:params.numFolds
        testIdx  = [edges(i):edges(i+1)-1 edges(i+params.numFolds):edges(i+params.numFolds+1)-1 ...
            edges(i+2*params.numFolds):edges(i+2*params.numFolds+1)-1 edges(i+3*params.numFolds):edges(i+3*params.numFolds+1)-1 ...
            edges(i+4*params.numFolds):edges(i+4*params.numFolds+1)-1];
        
        trainIdx = setdiff(1:size(deconvSignal,2),testIdx);
        
        d_data.iter{i}.deconvTrain = deconvSignal(:,trainIdx);
        d_data.iter{i}.deconvTest  = deconvSignal(:, testIdx);
        
        
        %% find cells without signal and delete
        d_data.iter{i}.cellIdx = params.cellIdx;
        d_data.iter{i} = findEmptyROIs(d_data.iter{i});
        
        [d_data.iter{i}, deleteIdxTrain, deleteIdxTest] = findEmptyIdx(d_data.iter{i});
        
        trainIdx(deleteIdxTrain) = [];
        testIdx(deleteIdxTest)  = [];
        
        %% calculate encoding model
        d_data.iter{i}.firingRate           = poissonGLMPos(d_data.iter{i}.deconvTrain, posBinned(trainIdx), params.smooth);
        
        % calculate liielihood: P(signal|position)
        likelihood_test     = llh_signalGivPos(d_data.iter{i}.deconvTest, d_data.iter{i}.firingRate);
        likelihood_train   = llh_signalGivPos(d_data.iter{i}.deconvTrain, d_data.iter{i}.firingRate);
        
        
        %% calculate posterior: P(position| signal) = P(signal|position)*P(position)
        % calculate prior P(position) --> flat prior, assume every position is equally likely
        prior        = ones(1,size(likelihood_train,2)).*1/size(likelihood_train,2);
        
        posterior_trainNotNorm = zeros(size(likelihood_train));
        posterior_train = zeros(size(likelihood_train));
        for t = 1:size(likelihood_train,1)
            posterior_trainNotNorm(t,:)             = likelihood_train(t,:).*prior;
            posterior_train(t,:)                    = posterior_trainNotNorm(t,:)./nansum(posterior_trainNotNorm(t,:));
            
            [~, maxInd]                             = nanmax(posterior_train(t,:));
            d_data.iter{i}.predPos_train(t)    = params.binGaps(maxInd);
        end
        
        % test
        posterior_testNotNorm = zeros(size(likelihood_test));
        posterior_test = zeros(size(likelihood_test));
        for t = 1:size(likelihood_test,1)
            posterior_testNotNorm(t,:)             = likelihood_test(t,:).*prior;
            posterior_test(t,:)                    = posterior_testNotNorm(t,:)./nansum(posterior_testNotNorm(t,:));
            
            [~, maxInd]                            = nanmax(posterior_test(t,:));
            d_data.iter{i}.predPos_test(t)         = params.binGaps(maxInd);
        end
        
        % go back from bin to actual position value
        d_data.iter{i}.realPos_train = params.binGaps(posBinned(trainIdx));
        d_data.iter{i}.realPos_test  = params.binGaps(posBinned(testIdx));
        
        d_data.rmse_test(i)  = plotDecodingPos(d_data.iter{i}.predPos_test, d_data.iter{i}.realPos_test, 'closeNotSave', [],[], max(params.binGaps));
        d_data.rmse_train(i) = plotDecodingPos(d_data.iter{i}.predPos_train, d_data.iter{i}.realPos_train, 'closeNotSave', [],[],max(params.binGaps));
        
    end
    
    % choose fold with median performance in test data:
    d_data.medianRMSE_test  = nanmedian(d_data.rmse_test(:));
    [~, k]                  = min(abs(d_data.rmse_test(:)-median(d_data.rmse_test(:))));
    
    % save corresponding train rmse's
    d_data.medianRMSE_train = d_data.rmse_train(k);
    
    d_data.meanRMSE_test    = nanmean(d_data.rmse_test);
    d_data.meanRMSE_train   = nanmean(d_data.rmse_train);
    
    d_data.medianTestFold   = k;
    d_data.f                = d_data.iter{k}.firingRate;
end

function sn = GetSn(Y, range_ff, method)
%% Estimate noise standard deviation

%% inputs:
%   Y: N X T matrix, fluorescence trace
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5];
end
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp';
end
if any(size(Y)==1)
    Y = reshape(Y, [], 1);
else
    Y = Y';
end

%% estimate the noise
[psdx, ff] = pwelch(Y, [],[],[], 1);
indf = and(ff>=range_ff(1), ff<=range_ff(2));
switch method
    case 'mean'
        sn=sqrt(mean(psdx(indf, :)/2));
    case 'median'
        sn=sqrt(median(psdx(indf,:)/2));
    case 'logmexp'
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));    
    otherwise
        fprintf('wrong method! use logmexp instead.\n'); 
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));
end
sn = sn';
end

function roiStat = getStats(sData, CH)

    if nargin < 2
        CH = 2; % Most of us are imaging on ch2?
    end
    
    dff = double(squeeze(sData.imdata.roiSignals(CH).dff));
    [nRois, nSamples] = size(dff);
    

    % Get max DFF of all rois.
    peakDff = max(dff, [], 2);
    
    
    % Get SNR of all Rois.
    signalToNoise = zeros(nRois, 1);
    noiseLevel = zeros(nRois, 1);
    for i = 1:nRois
        noiseLevel(i) = real(GetSn(dff(i, :)));
        signalToNoise(i) = snr(dff(i, :), ones(1,nSamples) * noiseLevel(i));
    end
    
    dffparams.smooth = params.smoothdata(dff, 2, 'movmean', 7); % params.smooth again with movmean
    
    isHigh = dffparams.smooth > noiseLevel;
    
    activityLevel = sum(isHigh, 2) ./ nSamples;
        
    roiStat = struct('peakDff', peakDff, ...
                     'signalToNoise', signalToNoise, ....
                     'activityLevel', activityLevel);
               
        
end


function [posBinned, signalBinned,params] = binData(pos, signal, params)

    timeBins = 1:params.timebinSize:length(pos);
    
    posBinned = zeros(length(timeBins)-1, 1);
    signalBinned = zeros(size(signal, 1), length(timeBins)-1);
    
    % this is redundant now but we tried differernt time binniings, and found
    % the best decoding with time bin = 1
    for t = 1:length(timeBins)-1
        posBinned(t) = nanmean(pos(timeBins(t):timeBins(t+1)-1));
        signalBinned(:, t) = nanmean(signal(:, timeBins(t):timeBins(t+1)-1), 2);
    end
    
    % remove all bins where animal is in position at beginniing oof the
    % track
    signalBinned(:, posBinned < params.deletePosThreshold) = [];
    posBinned(posBinned < params.deletePosThreshold) = [];
    
    params.binGaps(params.binGaps < params.deletePosThreshold-params.binSize/2) = [];

    for t = 1:length(posBinned)
        [~,ind] = min(abs(params.binGaps-posBinned(t)));
        posBinned(t) = ind;
    end    
    
end


function saveDecodingData(signalBinned, pos,params, dirct, fileName)

    d_data = main_decoding(signalBinned, pos, params);

    % Prepare for saving

    % Calculate confusion matrix
    classAccurTest = confusionMatrix(d_data.iter{d_data.medianTestFold}.predPos_test/params.binSize, d_data.iter{d_data.medianTestFold}.realPos_test/params.binSize, dirct, strcat('confMatr_test',fileName(8:end-4)), params.binGaps/params.binSize);
    classAccurTrain = confusionMatrix(d_data.iter{d_data.medianTestFold}.predPos_train/params.binSize, d_data.iter{d_data.medianTestFold}.realPos_train/params.binSize, dirct, strcat('confMatr_train',fileName(8:end-4)), params.binGaps/params.binSize);

    % Save the decoding data
    save(fullfile(dirct, fileName), 'd_data', '-v7.3');
end