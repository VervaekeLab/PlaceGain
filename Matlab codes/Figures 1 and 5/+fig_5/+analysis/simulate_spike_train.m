
% amplitude  deconv: 0.0004 max 0.06 use these in Amp values in code [0.01 0.05:0.05::0.5]
% with mean: amplitude 0.005 --> Amp = 0.05

% field width  deconv: 7cm max 20 cm use these in values  2 = 5 cm; 3 = 8 cm; 4 = 10 cm; 5 = 12 cm; 6 = 16 cm; 8 = 20 cm; 10 = 23;12 = 28 cm;
% with ca mean field width 10 cm --> sigma = 4 ;


sim.A             = 1.5;  % single AP DFF in % 
sim.Asigma        = [];   % spread of A around mean
sim.tauDecay      = 0.51; % indicator decay time in s / vectorizable for network of neurons
sim.tausigma      = [];   % spread of tau1 around mean
sim.tauRise       = 0.13; % 0.01; % onset/ rise time in s
sim.samplingRate  = 1000; % rate for calculation of calcium concentration and DFF in Hz
sim.frameRate     = 31;   %  % simulate imaging of DFF with this rate in Hz

sim.SNR           = 4;
sim.snr           = 4;
Amp               = [0.01 0.05:0.05:0.5]';
sigma             = 4 ;
offset            = [0 0.003,0.005:0.001:0.01,0.015 0.02];

[m,n] = ndgrid(Amp,offset);
Z = [m(:) n(:)];

% we load an experimental sData to get the position profile 
template_sData      = load(fullfile('...','/sData.mat'));
exp_sData           =  template_sData.sData;

% prepare saving directory
savedirct = fullfile('...','simulation');
if isempty(savedirct); mkdir(savedirct);end
    
downSampleDeltaD        =   1/10000; % Sampling rate in s based on the microscope frame clock.

for k =1:length(Z)
    
    sData{k}.behavior.wheelPos   = exp_sData.behavior.wheelPos(1:end);
    sData{k}.behavior.runSpeed   = exp_sData.behavior.runSpeed(1:end);
    sData{k}.behavior.runSpeedDs = exp_sData.behavior.runSpeedDs(1:end);
    position     = exp_sData.behavior.wheelPos(1:end);
    duration     = length(position)/1000; % simulation duration in s


    pos         = zeros(length(position),1);
    binGapsPos  = 0:max(position); 
    yCenter     = 1:length(binGapsPos);

    for t = 1:length(position)
        [~,ind] = min(abs(binGapsPos-position(t)));
        pos(t) = ind;
    end

    trialStart         =  find((diff(exp_sData.daqdata.waterValve) == 1)); %
    trialLength        = zeros(1,length(exp_sData.behavior.wheelPos)); % Initializes variable.
    for f = 1:numel(trialStart)-1
        trialLength(trialStart(f):trialStart(f+1))=f; % Assigns value i to sData.trials.trialLength corresponding to the trial number the index belongs to.
    end

    numNeurons = 30;

    x = 1:158;

    spikeRate = Z(k,2)+Z(k,1)*normpdf(x,78,sigma);
    
    for n = 1:numNeurons

            spikeRate_Temp = circshift(spikeRate,70+n*5);
            %     poisson spike train consistent across trials
            spkTimes = [];
            for i = 1:158
                idx                = find(pos == i);
                spkTimesTemp       = PoissonSpikeTrain(spikeRate_Temp(i),length(idx));
                spkTimes = [spkTimes; idx(spkTimesTemp)];
            end

        [dff, denoised, spikes, sData{k}]             = simulateSingleCa(sim, sData{k},exp_sData,unique(spkTimes),duration,position,trialLength);
        sData{k}.imdata.roiSignals(2).dff(n,:)        = dff;
        sData{k}.imdata.roiSignals(2).denoised(n,:)   = denoised;
        sData{k}.imdata.roiSignals(2).spikes(n,:)     = spikes;

    end

    savedData(savedirct,k, sData{k})
end


function savedData(savedirct,k, sData)
save(fullfile(savedirct, strcat('sData_',num2str(k))), '-struct', 'sData');

end

function [dff, denoised, spikes,sData] = simulateSingleCa(sim, sData,exp_sData, spkTimes,dur,position,trialLength)

x           = 0:1/sim.samplingRate:dur; x = x(1:end-1);
spkT        = x(spkTimes);
A           = sim.A;

% expected peak amplitude based on exponential rise and decay
PeakA       = A .* (sim.tauDecay./sim.tauRise.*(sim.tauDecay./sim.tauRise+1).^-(sim.tauRise./sim.tauDecay+1));
sdnoise     = PeakA./sim.snr; % determines how much Gaussian noise is added to the trace


% convolution of decay over all spikes (faster, same model calcium transient)
modelTransient      = spkTimes2Calcium(0,sim.tauRise,sim.A,sim.tauDecay,sim.samplingRate,dur);

spkVector           = zeros(1,numel(x));
spkVector(spkTimes) = 1;

% this is to simulate calcium signals
denoisedConv             = conv(spkVector,modelTransient);
denoised                 = denoisedConv(1:length(x));

% add noise
noisyDFF        = zeros(1,size(denoised,2));
whiteNoise      = sdnoise.*randn(1,size(denoised,2));
noisyDFF        = denoised + whiteNoise;

%% noisy DFF at frame rate
lowResT         = 1/sim.frameRate:1/sim.frameRate:max(x);
% %lowResT         = (1/S.frameRate:1/S.frameRate:max(x)) - (0.5 .* 1/S.frameRate);
 [~,idxList]     = findClosest(lowResT,x);
% 
%TimeDs          = x(idxList);
noisyDFFlowRes  = noisyDFF(idxList);

denoised        = denoised(idxList);
spikes          = spkVector;
dff             = noisyDFFlowRes;

sData.behavior.wheelPosDs = sData.behavior.wheelPos(idxList);
sData.stats.wheelCircum   = 158;

sData.behavior.wheelPosDs = position(idxList);%exp_sData.behavior.wheelPos(idxList);
sData.behavior.wheelPos   = position;
sData.behavior.runSpeedDs = exp_sData.behavior.runSpeed(idxList);
sData.behavior.runSpeedDs = exp_sData.behavior.runSpeed;
sData.trials.trialLength  = trialLength;
sData.trials.trialLengthDs  = trialLength(idxList);


end


function [y, x] = spkTimes2Calcium(spkT,tauRise,amp,tauDecay,frameRate,duration)

x = 1:(1/frameRate):duration;
y = (1-(exp(-(x-spkT)./tauRise))).*(amp*exp(-(x-spkT)./tauDecay));
y(x < spkT) = 0;
y(isnan(y)) = 0;

end

function spikeTimes = PoissonSpikeTrain(rate, dur)
% Generate Poisson Spike Train with firing rate and duration


    spikeTimes=[];
   
    % Generating spikes from a exponential distribution
    for t=1:dur
        c = rand;
        if  rate > c
            spikeTimes(end+1)   =   t;
        end
    end

end


function [d,ib] = findClosest(a,b)
% for each element of vector a find the nearest element in vector b based
% on absolute difference
% a and b must be row vectors
% d ... difference vector
% id ... indices of the differences recorded in d

% very efficiently implemented solution for this problem by Roger Stafford 
% (found on Matlab Central)

m = size(a,2); n = size(b,2);
[~,p] = sort([a,b]);
q = 1:m+n; q(p) = q;
t = cumsum(p>m);
r = 1:n; r(t(q(m+1:m+n))) = r;
s = t(q(1:m));
id = r(max(s,1));
iu = r(min(s+1,n));
[d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
ib = id+(it-1).*(iu-id);
end


