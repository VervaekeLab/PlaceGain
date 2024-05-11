dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/simulation/Alexander_and_Nitz/fig1c1d.mat';

load(dirct)

for i = 1:size(hpcfull,1)
    % Find peaks in the original activity matrix
    peak_activity_hpc(i)     =   max( hpcfull(i,:));
    baseline_activity_hpc(i ) =   prctile(hpcfull(i,:),5);

end


for i = 1:size(rscfull,1)
    % Find peaks in the original activity matrix
    peak_activity_rsc(i) =   max( rscfull(i,:));
    baseline_activity_rsc(i ) = prctile(rscfull(i,:),5);

end


fprintf('RSC peak activity: %6.2f pm %6.2f Hz \n',nanmean(peak_activity_rsc),nanstd(peak_activity_rsc));
fprintf('HPC peak activity: %6.2f pm %6.2f Hz \n',nanmean(peak_activity_hpc),nanstd(peak_activity_hpc));

fprintf('RSC background activity: %6.2f pm %6.5f Hz \n',nanmean(baseline_activity_rsc),nanstd(baseline_activity_rsc));
fprintf('HPC background  activity: %6.5f pm %6.5f Hz \n',nanmean(baseline_activity_hpc),nanstd(baseline_activity_hpc));


figure()
subplot(2,2,1)
h =histogram(peak_activity_hpc,'FaceColor',[0.5 0.5 0.5]);
xlabel('Peak activity (HPC)' )
ylabel('Cell #')
hold on
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;

box off 
subplot(2,2,2)
h =histogram(baseline_activity_hpc,'FaceColor',[0.5 0.5 0.5]);
xlabel('Baseline activity (HPC)' )
ylabel('Cell #')
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;

box off 
subplot(2,2,3)
h = histogram(peak_activity_rsc,'FaceColor',[0.5 0.5 0.5]);
xlabel('Peak activity (RSC)' )
ylabel('Cell #')
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;

box off 
subplot(2,2,4)
h =histogram(baseline_activity_rsc,'FaceColor',[0.5 0.5 0.5]);
xlabel('Baseline activity (RSC)' )
ylabel('Cell #')
box off
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;



