function [spatialSmoothness, bad_ic] = estimate_spatialsmoothnes(EEG, threshold_spatial)
% ESTIMATE_SPATIALSMOOTHNES Detects channel-like independent components
% based on spatial energy focalness (% of total variance in the peak electrode).
%
% Default threshold: 0.35 (flags ICs where >35% of total scalp power is in 1 channel)
%
% Variance:  < 3
% Kurtosis:  < 20-35
% Laplacian:  2000-3000

% Metric: Spatial Focalness (Energy Sparsity Ratio)
% Range is [1/n_channels, 1.0]. A pure brain dipole is typically < 0.15.
power_per_chan = EEG.icawinv.^2;
total_power    = sum(power_per_chan, 1);
max_power      = max(power_per_chan, [], 1);

spatialSmoothness = max_power ./ total_power;

% Flag ICs whose energy is disproportionately concentrated in one channel
bad_ic = find(spatialSmoothness > threshold_spatial);

% % Kurtosis
% spatialSmoothness = kurtosis(EEG.icawinv, 1, 1);
% bad_ic = find(spatialSmoothness > threshold_spatial);

end