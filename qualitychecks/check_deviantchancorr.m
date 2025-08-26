function [deviant_indices, fh] = check_deviantchancorr(CorrelationMatrices, subjects, options)
%FIND_DEVIANT_CORRELATION_MATRICES Identifies and visualizes outlier correlation matrices.
%
%   This function calculates the group average correlation matrix and identifies
%   individual matrices that deviate significantly from this average.
%
%   [deviant_indices, fh] = FIND_DEVIANT_CORRELATION_MATRICES(CorrelationMatrices, subjects, options)
%
%   INPUTS:
%   - CorrelationMatrices: An N-by-N-by-M matrix, where N is the number of
%                          channels and M is the number of subjects/trials.
%   - subjects:            A cell array of strings with M elements, containing
%                          the name or ID for each subject/trial.
%   - options:             (Optional) A struct with the following fields:
%       .metric:           The distance metric to use. Can be 'frobenius'
%                          (default) or 'mahalanobis' (not recommended!).
%       .z_threshold:      The Z-score threshold to define an outlier.
%                          (Default: 2)
%       .variance_to_keep: The percentage of variance to keep during PCA
%                          (only used for 'mahalanobis' metric).
%                          (Default: 0.95, i.e., 95%)
%
%   OUTPUTS:
%   - deviant_indices: A vector of indices for the identified outlier matrices.
%   - fh:              The handle to the generated figure.

% --- 1. Argument Parsing and Setup ---
if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'metric')
    options.metric = 'frobenius';
end
if ~isfield(options, 'z_threshold')
    options.z_threshold = 2;
end
if ~isfield(options, 'variance_to_keep')
    options.variance_to_keep = 0.95;
end

[n_channels, ~, n_subjects] = size(CorrelationMatrices);
if n_subjects ~= numel(subjects)
    error('The number of matrices must match the number of subject names.');
end

% --- 2. Calculate Group Average ---
CorrelationMatricesMean = trimmean(CorrelationMatrices,20,'round',3);

% --- 3. Calculate Deviation Score using the Chosen Metric ---
fprintf('Using ''%s'' distance metric to calculate deviation.\n', options.metric);

switch lower(options.metric)
    case 'mahalanobis'
        % --- PCA + Mahalanobis Distance Logic ---
        % First, vectorize the upper triangle of each correlation matrix.
        warning('Not a recommended method!');
        mask_upper_tri = triu(true(n_channels), 1);
        num_vars = sum(mask_upper_tri, 'all');
        vectorized_data = zeros(n_subjects, num_vars);
        for i = 1:n_subjects
            current_matrix = CorrelationMatrices(:,:,i);
            vectorized_data(i, :) = current_matrix(mask_upper_tri)';
        end

        % Check if dimensionality reduction is needed
        if n_subjects <= num_vars
            fprintf('Number of subjects (%d) is not greater than number of variables (%d).\n', n_subjects, num_vars);
            fprintf('Using PCA for dimensionality reduction.\n');
            [~, score, latent] = pca(vectorized_data);
            explained_variance = cumsum(latent) / sum(latent);
            num_components_to_keep = find(explained_variance >= options.variance_to_keep, 1, 'first');
            fprintf('Keeping %d components to explain %.2f%% of the variance.\n', num_components_to_keep, options.variance_to_keep * 100);
            data_for_mahal = score(:, 1:num_components_to_keep);
        else
            data_for_mahal = vectorized_data;
        end

        deviation_scores = mahal(data_for_mahal, data_for_mahal);
        xlabel_text = 'Deviation Score (Z-scored Mahalanobis)';

    case 'frobenius'
        % --- Frobenius Norm Logic ---
        deviation_scores = zeros(1, n_subjects);
        for i = 1:n_subjects
            diff_matrix = CorrelationMatrices(:,:,i) - CorrelationMatricesMean;
            deviation_scores(i) = norm(diff_matrix, 'fro');
        end
        xlabel_text = 'Deviation Score (Z-scored Frobenius Norm)';

    otherwise
        error("Invalid metric specified. Choose 'mahalanobis' or 'frobenius'.");
end

% --- 4. Identify Outliers ---
deviation_zscores = zscore(deviation_scores);
deviant_indices = find(abs(deviation_zscores) > options.z_threshold);

% Sort the outliers by how much they deviate, from most to least.
[~, sort_order] = sort(abs(deviation_zscores(deviant_indices)), 'descend');
deviant_indices = deviant_indices(sort_order);

fprintf('Found %d potential outliers with a Z-score > %.2f.\n', numel(deviant_indices), options.z_threshold);

% --- 5. Visualization ---
n_plots = numel(deviant_indices);
Ncol = 4;
% Nrow = ceil(n_plots / Ncol);
Nrow = n_plots;
Nrow = max(1,Nrow);

fh = figure('Color', 'w', 'Name', 'Correlation Matrix Outlier Check');
th = tiledlayout(1+Nrow,Ncol);
th.TileSpacing = 'compact'; th.Padding = 'compact';
title(th,{['Participants with larger covarance deviations, N = ' num2str(n_plots)],'Red flag: SVD1 is not brain-like and has >0.2'});

clim = [-max(abs(CorrelationMatricesMean(:))), max(abs(CorrelationMatricesMean(:)))];

% Plot the group average
nexttile(1);
imagesc(CorrelationMatricesMean, clim);
axis square; colorbar;
title('Group Average', 'FontWeight', 'bold');

% Plot the histogram of deviation scores
nexttile([1 3]);
histogram(deviation_zscores,15);
xlabel(xlabel_text); ylabel('Count');
title('Distribution of Deviation Scores', 'FontWeight', 'bold');
xlim([-4 4]);

% Plot the most deviant matrices
for i = 1:numel(deviant_indices)
    idx = deviant_indices(i);
    nexttile;
    imagesc(CorrelationMatrices(:,:,idx), clim);
    % imagesc(CorrelationMatrices(:,:,idx) - CorrelationMatricesMean, clim);
    axis square; colorbar;
    title(sprintf('%s\n(Z = %.2f)', subjects{idx}, deviation_zscores(idx)), 'Interpreter', 'none');

    tmp = CorrelationMatrices(:,:,idx) - CorrelationMatricesMean;
    [u, s, ~] = svd(tmp(1:128,1:128));
    s = diag(s);
    s = s / sum(s);
    for j = 1:3
        mytopoplot(u(:,j),[],['SVD' num2str(j) ': ' num2str(round(s(j),2))], nexttile); axis tight;
    end
end

end
