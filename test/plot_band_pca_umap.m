function plot_band_pca_umap(base_paths, output_dir, target_band)
% =========================================================================
% Script for PCA and UMAP clustering of targeted EEG frequency bands
% Dynamic cohort/visit colour gradient (Controls: Greens, Patients: Reds)
% Relative power normalization (channel-wise)
% =========================================================================
% INPUTS:
%   base_paths  - Cell array of strings to cohort directories
%                 e.g. {'F:\...\RS\ALS\T1', 'F:\...\RS\ALS\T2', 'F:\...\RS\CONTROL\T1'}
%   output_dir  - Output path for dashboard figures
%   target_band - 1x2 array defining frequency band in Hz (default [15 30])
% =========================================================================

if nargin < 3 || isempty(target_band)
    target_band = [15, 30]; % Default to Beta band
end

fprintf('\n==================================================================\n');
fprintf('Targeted Spectral Clustering: %.1f - %.1f Hz\n', target_band(1), target_band(2));
fprintf('==================================================================\n');

if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end

% 1. Identify all subject folders across all provided base paths
all_subjects  = {};
all_paths     = {};
cohort_groups = {}; % 'CONTROL' or 'ALS'
visit_nums    = []; % 1, 2, 3...
full_labels   = {}; % 'CONTROL_T1', 'ALS_T5', etc.

for p = 1:length(base_paths)
    current_path = regexprep(base_paths{p}, '[\\/]+$', '');
    [parent_dir, visit_str] = fileparts(current_path);
    [~, group_str] = fileparts(parent_dir);

    % Parse numerical visit number (e.g. 'T1' -> 1)
    v_tokens = regexp(visit_str, '\d+', 'match');
    if ~isempty(v_tokens)
        v_num = str2double(v_tokens{1});
    else
        v_num = 1;
    end

    c_name = sprintf('%s_%s', upper(group_str), upper(visit_str));

    d = dir(current_path);
    isub = [d(:).isdir];
    subj_folders = {d(isub).name}';
    subj_folders(ismember(subj_folders, {'.', '..', 'reports'})) = [];

    for s = 1:length(subj_folders)
        all_subjects{end+1}  = subj_folders{s};
        all_paths{end+1}     = current_path;
        cohort_groups{end+1} = upper(group_str);
        visit_nums(end+1)    = v_num;
        full_labels{end+1}   = c_name;
    end
end

num_subject = length(all_subjects);
band_features = cell(num_subject, 1);
valid_labels  = cell(num_subject, 1);
valid_groups  = cell(num_subject, 1);
valid_visits  = [];
valid_sub_ids = cell(num_subject, 1);
valid_mask    = false(num_subject, 1);

fprintf('Extracting relative power for %d datasets... \n', num_subject);

for i_subject = 1:num_subject
    sub_id       = all_subjects{i_subject};
    current_path = all_paths{i_subject};

    power_files = dir(fullfile(current_path, sub_id, '**', '*_rawpower.mat'));
    if isempty(power_files), continue; end

    file_path = fullfile(power_files(1).folder, power_files(1).name);

    try
        % Load into a temporary struct to prevent variable leakage across iterations
        tmp = load(file_path);

        if isfield(tmp, 'data_power')
            freqs = tmp.data_power.freq_post;
            psd_post_sub = tmp.data_power.psd_post(1:128, :);
        elseif isfield(tmp, 'freq_post') && isfield(tmp, 'psd_post')
            freqs = tmp.freq_post;
            psd_post_sub = tmp.psd_post(1:128, :);
        else
            continue;
        end

        % Identify frequency bins within the target band
        f_idx = freqs >= target_band(1) & freqs <= target_band(2);
        f_all = freqs < 49;

        % Extract 128 channels and calculate relative power per channel
        psd_band = psd_post_sub(:, f_idx);
        total_power = sum(psd_post_sub(:, f_all), 2);

        % Avoid division by zero
        total_power(total_power == 0) = eps;
        psd_rel = psd_band ./ total_power;

        % Flatten the 128 x Freqs matrix into a 1D feature vector
        band_flat = psd_rel(:)';

        band_features{i_subject} = band_flat;
        valid_labels{i_subject}  = full_labels{i_subject};
        valid_groups{i_subject}  = cohort_groups{i_subject};
        valid_visits(end+1)      = visit_nums(i_subject);
        valid_sub_ids{i_subject} = sprintf('%s_%s', full_labels{i_subject}, sub_id);
        valid_mask(i_subject)    = true;

    catch ME
        warning('Failed processing %s: %s', sub_id, ME.message);
    end
end

% Filter to valid datasets
band_features   = band_features(valid_mask);
valid_labels    = valid_labels(valid_mask);
valid_groups    = valid_groups(valid_mask);
valid_sub_ids   = valid_sub_ids(valid_mask);
valid_sub_count = length(band_features);

if valid_sub_count < 3
    error('Not enough valid subjects to run clustering.');
end

X_features = cell2mat(band_features);
X_features(isnan(X_features)) = nanmedian(X_features(:));

fprintf('Retained %d datasets. Feature length: %d per subject.\n', ...
    valid_sub_count, size(X_features, 2));

% -------------------------------------------------------------------------
% 1. Dimensionality Reduction (PCA)
% -------------------------------------------------------------------------
X_norm = (X_features - mean(X_features, 1)) ./ std(X_features, 0, 1);

fprintf('Running PCA...\n');
[~, score, ~, ~, explained] = pca(X_norm);

cum_var = cumsum(explained);
num_pcs_to_keep = find(cum_var >= 90, 1);
if isempty(num_pcs_to_keep) || num_pcs_to_keep > 30
    num_pcs_to_keep = min(30, size(score, 2));
end

% -------------------------------------------------------------------------
% 2. PCA-Guided UMAP
% -------------------------------------------------------------------------
n_neighbors = min(75, floor(valid_sub_count / 3));
fprintf('Running UMAP on the top %d Principal Components...\n', num_pcs_to_keep);

[Y_umap, ~] = run_umap(score(:, 1:num_pcs_to_keep), ...
    'n_neighbors', n_neighbors, ...
    'min_dist', 0.1, ...
    'metric', 'euclidean', ...
    'method', 'MATLAB', ...       % Native MATLAB gradient descent (avoids MEX checks)
    'ask_to_save', false, ...     % Suppress dialogs asking to save templates
    'verbose', 'text');           % Text console feedback without Java GUI frames

% -------------------------------------------------------------------------
% 3. Dynamic Colour Palette Assignment
% -------------------------------------------------------------------------
% Controls: Light green -> Dark green across follow-ups
% Patients: Light red   -> Dark red   across follow-ups

unique_labels = unique(valid_labels);
num_categories = length(unique_labels);
colour_map = zeros(num_categories, 3);

ctrl_visits = sort(unique(valid_visits(contains(valid_groups, 'CONTROL'))));
als_visits  = sort(unique(valid_visits(~contains(valid_groups, 'CONTROL'))));

for k = 1:num_categories
    lbl = unique_labels{k};
    v_tok = regexp(lbl, '\d+', 'match');
    v_val = str2double(v_tok{1});

    if contains(lbl, 'CONTROL')
        v_idx = find(ctrl_visits == v_val, 1);
        n_v = max(length(ctrl_visits), 1);
        if n_v == 1, frac = 0.4; else, frac = (v_idx - 1) / (n_v - 1); end

        c_light = [0.65, 0.90, 0.65];
        c_dark  = [0.05, 0.40, 0.15];
        colour_map(k, :) = (1 - frac) * c_light + frac * c_dark;
    else
        v_idx = find(als_visits == v_val, 1);
        n_v = max(length(als_visits), 1);
        if n_v == 1, frac = 0.4; else, frac = (v_idx - 1) / (n_v - 1); end

        c_light = [0.95, 0.60, 0.60];
        c_dark  = [0.50, 0.05, 0.10];
        colour_map(k, :) = (1 - frac) * c_light + frac * c_dark;
    end
end

% -------------------------------------------------------------------------
% Visualisation 1: PCA Scatter
% -------------------------------------------------------------------------
fh1 = figure('Name', 'PCA of Target Band', 'Color', 'w', 'Position', [100, 200, 850, 650]);
hold on; grid on;

for k = 1:num_categories
    idx = strcmp(valid_labels, unique_labels{k});
    scatter(score(idx, 1), score(idx, 2), 65, colour_map(k, :), 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', unique_labels{k});
end

xlabel(sprintf('PC1 (%.1f%% Variance)', explained(1)), 'FontWeight', 'bold');
ylabel(sprintf('PC2 (%.1f%% Variance)', explained(2)), 'FontWeight', 'bold');
title(sprintf('PCA of Relative Power (%.1f-%.1f Hz, N=%d)', target_band(1), target_band(2), valid_sub_count), 'FontSize', 13);
legend('Location', 'best', 'Interpreter', 'none');

plotX = 22; plotY = 16;
set(fh1, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh1, fullfile(output_dir, sprintf('Dashboard_PCA_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 2: PCA-Guided UMAP
% -------------------------------------------------------------------------
fh2 = figure('Name', 'PCA-Guided UMAP', 'Color', 'w', 'Position', [980, 200, 850, 650]);
hold on; grid on;

for k = 1:num_categories
    idx = strcmp(valid_labels, unique_labels{k});
    scatter(Y_umap(idx, 1), Y_umap(idx, 2), 65, colour_map(k, :), 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', unique_labels{k});
end

xlabel('UMAP Dimension 1', 'FontWeight', 'bold');
ylabel('UMAP Dimension 2', 'FontWeight', 'bold');
title(sprintf('PCA-Guided UMAP (%.1f-%.1f Hz, Neighbors: %d)', target_band(1), target_band(2), n_neighbors), 'FontSize', 13);
legend('Location', 'best', 'Interpreter', 'none');

set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(output_dir, sprintf('Dashboard_UMAP_PCA_Guided_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

fprintf('Finished. Both dashboards saved to %s\n', output_dir);
end