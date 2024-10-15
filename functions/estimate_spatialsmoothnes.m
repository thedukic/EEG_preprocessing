function [spatialSmoothness, badIC, spatialTreshold] = estimate_spatialsmoothnes(EEG)

% Empirical tresholds
spatialTreshold = [3 3000];

% Calculate the spatial variance
theseICs = EEG.icawinv;
spatialSmoothness = NaN(2,size(EEG.icawinv,2));

% Spatial variance, low for localised ICs
thisICNorm = theseICs ./ max(abs(theseICs),[],1);
spatialSmoothness(1,:) = 100*var(thisICNorm);

% Laplacian variance, high for localised ICs
spatialLaplacian = estimate_laplacian(theseICs,EEG.chanlocs,1);
spatialSmoothness(2,:) = var(spatialLaplacian,0,1);

% figure; scatter(spatialSmoothness(1,:),spatialSmoothness(2,:));
% xlabel('Spatial variance'); ylabel('Laplacian variance');

% Tresholding
badIC = false(size(spatialSmoothness));
badIC(1,:) = spatialSmoothness(1,:) < spatialTreshold(1);
badIC(2,:) = spatialSmoothness(2,:) > spatialTreshold(2); % or 1000?
badIC = find(all(badIC,1));

end