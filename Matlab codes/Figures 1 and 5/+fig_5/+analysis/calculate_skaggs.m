
function [SI, specificity]  = calculate_skaggs(occupancy_map_2D, activity_map_2D)
%% CALC_SPATIALINFORMATION 
% Calculates spatial information according to the method used by Skaggs et al (1993)
%
% INPUT
%   occupancy_map: A matrix where each row corresponds to a lap of
%       running and each column is a binned position on the wheel. Each
%       element contains the number of samples the animal was in the bin.
%   activity_map_2D: A matrix where each row corresponds to a lap of
%       running and each column is a binned position on the wheel. Each
%       element contains summed activity of a signal (fex deconvolved) 
%       when the animal was in the bin.
%
% OUTPUT
%   SI: The spatial information as bits. 

%% Parameters

% bin_size = 1; %cm
n_bins = size(activity_map_2D,2); % The track number of bins

% Down project maps from 2D to 1D
occupancy_map_1D = sum(occupancy_map_2D,1);
% activity_map_1D = sum(activity_map_2D,1);

%% Calculate spatial information
% Create a position corrected activity map
% position_activity_map_1D = activity_map_1D./occupancy_map_1D;
mean_firing_of_bin_1D = nanmean(activity_map_2D);

% Smooth position activity map using a Gaussian window if deconv is used
% position_activity_map_1D = smoothdata(position_activity_map_1D,'gaussian',7); % position_activity_map_1D; %

% Create occupancy proability array
occupancy_probability = zeros(size(occupancy_map_1D));
for p = 1:length(occupancy_map_1D)
    occupancy_probability(p) = occupancy_map_1D(p)/nansum(occupancy_map_1D);
end

SI = NaN(n_bins,1);
for bin = 1:n_bins
    SI(bin) = occupancy_probability(bin) * (mean_firing_of_bin_1D(bin)/nanmean(mean_firing_of_bin_1D)) * ...
        log2(mean_firing_of_bin_1D(bin)/nanmean(mean_firing_of_bin_1D));
end


SI = nansum(SI)./0.14;
% information rate divided by overall mean firing rate leads to information
% per spike/specificity 
specificity  = SI/mean(mean_firing_of_bin_1D);



end