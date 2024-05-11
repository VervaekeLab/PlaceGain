dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/simulation/decoding';
dirct_sData ='/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/simulation/spikes';

% load simulated spikes
Amp               = [0.01 0.05:0.05:0.5]/0.001;
offset            = [0 0.003,0.005:0.001:0.01,0.015 0.02]/0.001;
sigma = 4;

[m,n] = ndgrid(Amp,offset);
Z = [m(:) n(:)];

x = 1:158;
l =1;
for k =1:length(Z)
    spikeRate = Z(k,2)+Z(k,1)*normpdf(x,78,sigma);
    spike_rate_vector(k,1) = max(spikeRate);
    spike_rate_vector(k,2) = min(spikeRate);
end

examples          =16;

for k = 1:length(examples)
     sData = load(strcat(dirct_sData,'/sData_',num2str(examples(k)),'.mat'));
     signal =sData.imdata.roiSignals(2).spikes(1:5:end,:);
end


position = sData.behavior.wheelPos;
speed =sData.behavior.runSpeed;


fig          = figure('Units','normalized','Position',[0 0 0.4 0.7],'color','w');

ax_signal_1     =  axes('Units','normalized','Position', [0.05 0.82 0.9 0.10]);
ax_signal_2     =  axes('Units','normalized','Position', [0.05 0.69 0.9 0.10]);
ax_signal_3     =  axes('Units','normalized','Position', [0.05 0.56 0.9 0.10]);
ax_signal_4     =  axes('Units','normalized','Position', [0.05 0.43 0.9 0.10]);
ax_signal_5     =  axes('Units','normalized','Position', [0.05 0.3 0.9 0.10]);

ax_behavior   =  axes('Units','normalized','Position', [0.05 0.05 0.9 0.2]);

stem(ax_signal_1,signal(2,1.0e+05 *1.5494:1.0e+05 *1.9312), 'Marker', 'none', 'Color', 'k', 'LineWidth', 0.1);
hold on
stem(ax_signal_2,signal(3,1.0e+05 *1.5494:1.0e+05 *1.9312), 'Marker', 'none', 'Color', 'k', 'LineWidth', 0.1);
stem(ax_signal_3,signal(4,1.0e+05 *1.5494:1.0e+05 *1.9312), 'Marker', 'none', 'Color', 'k', 'LineWidth', 0.1);
stem(ax_signal_4,signal(5,1.0e+05 *1.5494:1.0e+05 *1.9312), 'Marker', 'none', 'Color', 'k', 'LineWidth', 0.1);
stem(ax_signal_5,signal(6,1.0e+05 *1.5494:1.0e+05 *1.9312), 'Marker', 'none', 'Color', 'k', 'LineWidth', 0.1);


box(ax_signal_1, 'off'); axis(ax_signal_1, 'off');
box(ax_signal_2, 'off');axis(ax_signal_2, 'off');
box(ax_signal_3, 'off');axis(ax_signal_3, 'off');
box(ax_signal_4, 'off');axis(ax_signal_4, 'off');
box(ax_signal_5, 'off');axis(ax_signal_5, 'off');


plot(ax_behavior,position(1.0e+05 *1.5494:1.0e+05 *1.9312), 'LineWidth',1.5)

