dirct = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure5/analysis/simulation/decoding';


Amp               = [0.01 0.05:0.05:0.5]/0.001;
offset            = [0 0.003,0.005:0.001:0.01,0.015 0.02]/0.001;

[m,n] = ndgrid(Amp,offset);
Z = [m(:) n(:)];
sigma = 4;

x = 1:158;
for k =1:length(Amp)
    for i =1:length(offset)    
        spikeRate = offset(i)+Amp(k)*normpdf(x,78,sigma);
        X(k,i) = max(spikeRate);
        Y(k,i) = min(spikeRate);
    end   
end


for k = 1:length(Z)
    
    load(strcat(dirct,'/d_data_',num2str(k),'.mat'));
    rmse(k)   = nanmean(d_data.rmse_val);
end

l = 1;
for k = 1:length(offset)
    for i = 1:length(Amp)
        rmsegrid(i,k) = rmse(l);        
        l = l+1; 
    end
end


X_orig = X;
Y_orig = Y;
rmsegrid_orig = rmsegrid;

yq = offset ;
xq = 50*ones(1,length(yq));
vq = griddata(X,Y,rmsegrid, xq',yq');

Y(8:end,:)=[];
X(8:end,:)=[];
rmsegrid(8:end,:)=[];
rmsegrid(8,:)= vq';
X(8,:) = xq;
Y(8,:) = yq;

% X(:,1) = [];
% Y(:,1) = [];
% rmsegrid(:,1) = [];
% interpolation of missing value
rmsegrid(end,1) =rmsegrid(end-1,1)-(rmsegrid(end-2,1) -rmsegrid(end-1,1));
fig = figure();
a =surf(X,Y,rmsegrid,'FaceLighting','gouraud',...
    'MeshStyle','column',...
    'SpecularColorReflectance',0,...
    'SpecularExponent',5,...
    'SpecularStrength',0.2,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'FaceAlpha',0.2,...
    'FaceColor',[0.8 0.8 0.8],...
    'EdgeAlpha',0.2);%'EdgeColor',[0 0 0],'LineStyle','none','FaceAlpha',0.3);
colormap('gray')
a.CData = 36.1241*ones(size(a.CData));

ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
hold on
xlabel('Peak firing rate (Hz)')
zlabel('Decoding error (cm)')
ylabel('Out-of-field firing rate (Hz)')




for i = 0:5:50
    yq = 0.0001:20;
    xq = i*ones(1,size(yq,2));
    vq = griddata(X,Y,rmsegrid, xq',yq');
    plot3( xq',yq',vq,'Color',[0.5 0.5 0.5],'LineWidth',1)
end



% find_values_Andrew_Nitz;
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

baseline_activity_rsc = baseline_activity_rsc+0.0001;
plane_value_rsc= griddata(X_orig,Y_orig,rmsegrid_orig,peak_activity_rsc,baseline_activity_rsc);
scatter3(peak_activity_rsc,baseline_activity_rsc,plane_value_rsc,50,'filled','r','MarkerFaceAlpha',0.5);

baseline_activity_hpc = baseline_activity_hpc+0.0001;
plane_value_hpc= griddata(X_orig,Y_orig,rmsegrid_orig,peak_activity_hpc,baseline_activity_hpc);
scatter3(peak_activity_hpc,baseline_activity_hpc,plane_value_hpc,50,'filled','b','MarkerFaceAlpha',0.5);


% example RSC mean firing rate
figure()
subplot(2,1,1)
plot(rscfull(76,1:36),'r','LineWidth',1.52) % 76,83, 102,110

subplot(2,1,2)
plot(rscfull(141,1:36),'r','LineWidth',1.52) % 76,83, 102,110


figure()
plot(rscfull(53,1:36),'r','LineWidth',1.52) % 76,83, 102,110




