function plot_band_pca_tsne(base_paths, output_dir, target_band)
% =========================================================================
% Script for PCA and t-SNE clustering of targeted EEG frequency bands
% =========================================================================
% INPUTS:
%   base_paths  - Cell array of strings to cohort directories
%   output_dir  - String specifying where to save the dashboards
%   target_band - 1x2 array defining the frequency band in Hz (e.g. [13 30])
% =========================================================================

if nargin < 3
    target_band = [13, 30]; % Default to Beta band if not specified
end

fprintf('\n==================================================================\n');
fprintf('Targeted Spectral Clustering: %.1f - %.1f Hz\n', target_band(1), target_band(2));
fprintf('==================================================================\n');

if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end

% 1. Identify all subject folders
all_subjects = {};
all_paths    = {};
cohort_names = {};

for p = 1:length(base_paths)
    current_path = regexprep(base_paths{p}, '[\\/]+$', '');
    [parent_dir, visit_name] = fileparts(current_path);
    [~, cohort_main] = fileparts(parent_dir);
    c_name = sprintf('%s_%s', cohort_main, visit_name);

    d = dir(current_path);
    isub = [d(:).isdir];
    subj_folders = {d(isub).name}';
    subj_folders(ismember(subj_folders, {'.', '..', 'reports'})) = [];

    for s = 1:length(subj_folders)
        all_subjects{end+1} = subj_folders{s};
        all_paths{end+1}    = current_path;
        cohort_names{end+1} = c_name;
    end
end

num_subject = length(all_subjects);
band_features = cell(num_subject, 1);
valid_labels  = cell(num_subject, 1);
valid_sub_ids = cell(num_subject, 1);
valid_mask    = false(num_subject, 1);

fprintf('Extracting targeted power for %d datasets... \n', num_subject);

for i_subject = 1:num_subject
    sub_id       = all_subjects{i_subject};
    current_path = all_paths{i_subject};
    c_name       = cohort_names{i_subject};

    power_files = dir(fullfile(current_path, sub_id, '**', '*_rawpower.mat'));
    if isempty(power_files), continue; end

    file_path = fullfile(power_files(1).folder, power_files(1).name);

    try
        % 1. Load data into a temporary struct to prevent workspace leakage
        tmp = load(file_path);

        % 2. Check how the data was saved (inside a struct vs directly in the file)
        if isfield(tmp, 'data_power')
            freqs = tmp.data_power.freq_post;
            psd_post_sub = tmp.data_power.psd_post(1:128, :);
        elseif isfield(tmp, 'freq_post') && isfield(tmp, 'psd_post')
            freqs = tmp.freq_post;
            psd_post_sub = tmp.psd_post(1:128, :);
        else
            error('Neither ''data_power'' struct nor ''psd_post'' variable found.');
        end

        % 3. Identify frequency bins within the target band
        f_idx = freqs >= target_band(1) & freqs <= target_band(2);

        if sum(f_idx) == 0
            warning('No frequency bins found in the %.1f-%.1f Hz range for %s', ...
                target_band(1), target_band(2), sub_id);
            continue;
        end

        % Extract the 128 channels for the target band only
        psd_band = psd_post_sub(:, f_idx);

        % Log-transform to decibels
        psd_db = 10 * log10(psd_band);

        % Flatten the matrix into a 1D vector (Chan x Freqs)
        band_flat = psd_db(:)';

        band_features{i_subject} = band_flat;
        valid_labels{i_subject}  = c_name;
        valid_sub_ids{i_subject} = sprintf('%s_%s', c_name, sub_id);
        valid_mask(i_subject)    = true;

    catch ME
        warning('Failed to process data for %s: %s', sub_id, ME.message);
    end
end

% Filter to valid entries
band_features   = band_features(valid_mask);
valid_labels    = valid_labels(valid_mask);
valid_sub_ids   = valid_sub_ids(valid_mask);
valid_sub_count = length(band_features);

if valid_sub_count < 3
    error('Not enough valid subjects to run clustering.');
end

X_features = cell2mat(band_features);
X_features(isnan(X_features)) = nanmedian(X_features(:));

fprintf('Retained %d subjects. Feature length reduced to %d per subject.\n', ...
    valid_sub_count, size(X_features, 2));

% -------------------------------------------------------------------------
% 1. Principal Component Analysis (PCA)
% -------------------------------------------------------------------------
% Z-score standardisation before PCA is critical for PSD data
X_norm = (X_features - mean(X_features, 1)) ./ std(X_features, 0, 1);

fprintf('Running PCA...\n');
[coeff, score, latent, tsquared, explained] = pca(X_norm);

% Visualisation 1: Pure PCA Scatter
fh1 = figure('Name', 'PCA of Target Band', 'Color', 'w', 'Position', [100, 200, 800, 600]);
hold on; grid on;

unique_cohorts = unique(valid_labels);
colours = lines(length(unique_cohorts));

for c = 1:length(unique_cohorts)
    idx = strcmp(valid_labels, unique_cohorts{c});
    scatter(score(idx, 1), score(idx, 2), 70, colours(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8, 'DisplayName', unique_cohorts{c});
end

xlabel(sprintf('PC1 (%.1f%% Variance)', explained(1)), 'FontWeight', 'bold');
ylabel(sprintf('PC2 (%.1f%% Variance)', explained(2)), 'FontWeight', 'bold');
title(sprintf('PCA of Log-Power (%.1f-%.1f Hz)', target_band(1), target_band(2)), 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'none');

plotX = 22; plotY = 16;
set(fh1, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh1, fullfile(output_dir, sprintf('Dashboard_PCA_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% 2. PCA-Guided t-SNE
% -------------------------------------------------------------------------
% Determine how many PCs to feed into t-SNE (retain 90% variance or max 30)
cum_var = cumsum(explained);
num_pcs_to_keep = find(cum_var >= 90, 1);
if isempty(num_pcs_to_keep) || num_pcs_to_keep > 30
    num_pcs_to_keep = min(30, size(score, 2));
end

fprintf('Running t-SNE on the top %d Principal Components...\n', num_pcs_to_keep);

rng('default');
% With N=275, perplexity between 30 and 50 is optimal.
tsne_perplexity = min(10, max(5, floor(valid_sub_count / 4)));

% Feed the PCA scores (not the raw features) into t-SNE
Y_tsne = tsne(score(:, 1:num_pcs_to_keep), 'Standardize', true, 'Perplexity', tsne_perplexity);

% Visualisation 2: PCA-Guided t-SNE
fh2 = figure('Name', 'PCA-Guided t-SNE', 'Color', 'w', 'Position', [950, 200, 800, 600]);
hold on; grid on;

for c = 1:length(unique_cohorts)
    idx = strcmp(valid_labels, unique_cohorts{c});
    scatter(Y_tsne(idx, 1), Y_tsne(idx, 2), 70, colours(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8, 'DisplayName', unique_cohorts{c});
end

xlabel('t-SNE Dimension 1', 'FontWeight', 'bold');
ylabel('t-SNE Dimension 2', 'FontWeight', 'bold');
title(sprintf('PCA-Guided t-SNE (%.1f-%.1f Hz, Perplexity: %d)', ...
    target_band(1), target_band(2), tsne_perplexity), 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'none');

set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(output_dir, sprintf('Dashboard_tSNE_PCA_Guided_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

fprintf('Finished. Plots saved to %s\n', output_dir);
end