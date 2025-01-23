function [spatialSmoothness, badIC, spatialTreshold] = estimate_spatialsmoothnes(EEG)

% Empirical tresholds
spatialTreshold = [3 2000];

% Calculate the spatial variance
theseICs = EEG.icawinv;
spatialSmoothness = NaN(2,size(EEG.icawinv,2));

% Spatial variance, low for localised ICs
thisICNorm = theseICs ./ max(abs(theseICs),[],1);
spatialSmoothness(1,:) = 100 * var(thisICNorm,0,1);

% Laplacian variance, high for localised ICs
spatialLaplacian = estimate_laplacian(theseICs,EEG.chanlocs,1);
spatialSmoothness(2,:) = var(spatialLaplacian,0,1);

% Tresholding
badIC = false(size(spatialSmoothness));
badIC(1,:) = spatialSmoothness(1,:) < spatialTreshold(1);
badIC(2,:) = spatialSmoothness(2,:) > spatialTreshold(2);
badIC = find(all(badIC,1));

% figure; scatter(spatialSmoothness(1,:),spatialSmoothness(2,:));
% xlabel('Spatial variance'); ylabel('Laplacian variance');

% badIC = spatialSmoothness(1,:) < 3;
% badIC = spatialSmoothness(2,:) > 3000;
% badIC = find(badIC);

% figure;
% for i = 1:length(badIC)
%     mytopoplot(theseICs(:,badIC(i)), [],num2str(badIC(i)),nexttile);
% end

end