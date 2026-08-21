function [spatialSmoothness, badIC] = estimate_spatialsmoothnes(EEG, spatialTreshold)
% Variance:  < 3
% Kurtosis:  < 20-35
% Laplacian:  2000-3000

% Calculate the spatial variance
theseICs = EEG.icawinv;
% spatialSmoothness = NaN(2,size(EEG.icawinv,2));

% % Metric 1: Spatial Variance:
% % By normalising the weights so the peak is 1,
% % a genuine bad channel Independent Component will look like [1, 0.01, 0.02, 0, 0.01...].
% % Taking the variance of this array yields a very low number because
% % the vast majority of the values are near zero. Conversely, a real brain signal
% % spreads across the scalp, meaning more values are non-zero, resulting in a higher variance.
% thisICNorm = theseICs ./ max(abs(theseICs),[],1);
% % spatialSmoothness(1,:) = 100 * var(thisICNorm,0,1); % < 3
% spatialSmoothness(1, :) = kurtosis(thisICNorm, 1, 1); % > 35
%
% % a = 100 * var(thisICNorm,0,1);
% % b = kurtosis(thisICNorm, 1, 1);
% % figure; scatter(a, b);
% % xlabel('Spatial variance'); ylabel('Kurtosis');
%
% % Metric 2: Laplacian Variance: The spatial Laplacian calculates the local
% % curvature or sharpness between neighbouring electrodes. A single popping
% % channel creates an extreme spike on the scalp map. This sharp gradient
% % results in a massive Laplacian value at that specific coordinate and
% % near-zero values elsewhere, leading to a very high overall variance.
% spatialLaplacian = estimate_laplacian(theseICs,EEG.chanlocs,1);
% spatialSmoothness(2,:) = var(spatialLaplacian,0,1);
%
% % Tresholding
% badIC = false(size(spatialSmoothness));
% badIC(1,:) = spatialSmoothness(1,:) > spatialTreshold(1);
% badIC(2,:) = spatialSmoothness(2,:) > spatialTreshold(2);
% % badIC = find(all(badIC,1));
% % badIC = find(any(badIC,1));
% badIC = find(badIC(1,:));

% Kurtosis
thisICNorm = theseICs ./ max(abs(theseICs),[],1);
spatialSmoothness = kurtosis(thisICNorm, 1, 1);
badIC = find(spatialSmoothness > spatialTreshold);

% figure; scatter(spatialSmoothness(1,:),spatialSmoothness(2,:));
% xlabel('Spatial variance'); ylabel('Laplacian variance');

% figure;
% for i = 1:length(badIC)
%     mytopoplot(theseICs(:,badIC(i)), [],num2str(badIC(i)),nexttile);
% end

end