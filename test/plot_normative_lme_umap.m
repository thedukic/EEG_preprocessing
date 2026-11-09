function plot_normative_lme_umap(base_paths, output_dir, table_path, target_band)
% =========================================================================
% Continuous Spatio-Spectral NMF & Normative LME Pipeline
% 1. Retains all continuous frequency bins without sub-band averaging.
% 2. Performs continuous spatio-spectral NMF (Channels x Frequencies).
% 3. Fits normative LME models (REML) on Healthy Controls per component.
% 4. Generates side-by-side Topographies and Continuous Spectral Tuning curves.
% 5. Projects NMF normative residuals into Euclidean UMAP space.
% =========================================================================

if nargin < 4 || isempty(target_band)
    target_band = [2, 30]; % Broad physiological window (avoids arbitrary cuts)
end

if nargin < 3 || isempty(table_path)
    table_path = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-10-13.txt';
end

if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end

% -------------------------------------------------------------------------
% 1. Load Clinical Metadata
% -------------------------------------------------------------------------
fprintf('\nLoading Clinical Data from %s...\n', table_path);
Table1 = readtable(table_path, 'VariableNamingRule', 'preserve');
Table1.ID  = categorical(cellstr(Table1.ID));
Table1.SEX = categorical(cellstr(Table1.SEX));

% -------------------------------------------------------------------------
% 2. Parse Subject Directories
% -------------------------------------------------------------------------
all_subjects  = {};
all_paths     = {};
cohort_groups = {};
visit_nums    = [];
full_labels   = {};

for p = 1:length(base_paths)
    current_path = regexprep(base_paths{p}, '[\\/]+$', '');
    [parent_dir, visit_str] = fileparts(current_path);
    [~, group_str] = fileparts(parent_dir);

    v_tokens = regexp(visit_str, '\d+', 'match');
    v_num = str2double(v_tokens{1});
    if isempty(v_tokens), v_num = 1; end

    c_name = sprintf('%s_%s', upper(group_str), upper(visit_str));
    d = dir(current_path);
    subj_folders = {d([d(:).isdir]).name}';
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
raw_sub_ids   = cell(num_subject, 1);
valid_mask    = false(num_subject, 1);
sub_freqs     = [];

fprintf('Extracting continuous spectral power (%.1f - %.1f Hz) for %d datasets... \n', ...
    target_band(1), target_band(2), num_subject);

for i_subject = 1:num_subject
    sub_id       = all_subjects{i_subject};
    current_path = all_paths{i_subject};

    power_files = dir(fullfile(current_path, sub_id, '**', '*_rawpower.mat'));
    if isempty(power_files), continue; end

    try
        tmp = load(fullfile(power_files(1).folder, power_files(1).name));

        mask_channel = 1:128;
        if isfield(tmp, 'data_power')
            freqs = tmp.data_power.freq_post;
            psd_post = tmp.data_power.psd_post(mask_channel, :);
        elseif isfield(tmp, 'freq_post') && isfield(tmp, 'psd_post')
            freqs = tmp.freq_post;
            psd_post = tmp.psd_post(mask_channel, :);
        else
            continue;
        end

        f_idx = freqs >= target_band(1) & freqs <= target_band(2);
        f_all = freqs >= 1 & freqs <= 45; % Bounded denominator to eliminate drift/noise

        if isempty(sub_freqs)
            sub_freqs = freqs(f_idx);
        end


        psd_band = psd_post(:, f_idx);
        total_power = sum(psd_post(:, f_all), 2);
        total_power(total_power == 0) = eps;

        % Retain continuous frequency bins (128 channels x n_freqs) without averaging
        psd_rel = psd_band ./ total_power;
        band_flat = psd_rel(:)';

        band_features{i_subject} = band_flat;
        valid_labels{i_subject}  = full_labels{i_subject};
        valid_groups{i_subject}  = cohort_groups{i_subject};
        valid_visits(end+1)      = visit_nums(i_subject);
        raw_sub_ids{i_subject}   = sub_id;
        valid_mask(i_subject)    = true;

    catch ME
        warning('Failed processing %s: %s', sub_id, ME.message);
    end
end

band_features = band_features(valid_mask);
valid_labels  = valid_labels(valid_mask);
valid_groups  = valid_groups(valid_mask);
raw_sub_ids   = raw_sub_ids(valid_mask);
valid_sub_count = length(band_features);

% -------------------------------------------------------------------------
% 3. Match Clinical Metadata & Compute Longitudinal Variables
% -------------------------------------------------------------------------
Table1.DATE = datetime(Table1.("DATE: VISIT"));

matched_age_base = NaN(valid_sub_count, 1);
matched_sex      = categorical(repmat({''}, valid_sub_count, 1));
matched_hand     = categorical(repmat({''}, valid_sub_count, 1));
matched_time     = NaN(valid_sub_count, 1);
clinical_mask    = false(valid_sub_count, 1);

for i = 1:valid_sub_count
    subj_id = categorical(cellstr(raw_sub_ids{i}));
    current_visit = valid_visits(i);

    idx_current = find(Table1.ID == subj_id & Table1.VISIT == current_visit);
    idx_base    = find(Table1.ID == subj_id & Table1.VISIT == 1);

    if ~isempty(idx_current) && ~isempty(idx_base)
        matched_age_base(i) = Table1.AGE(idx_base(1));
        matched_sex(i)      = Table1.SEX(idx_base(1));

        if strcmpi(Table1.Handedness(idx_base(1)), 'NA') || isempty(Table1.Handedness(idx_base(1)))
            hand = 'R';
        else
            hand = Table1.Handedness(idx_base(1));
        end
        matched_hand(i) = hand;

        date_base    = Table1.DATE(idx_base(1));
        date_current = Table1.DATE(idx_current(1));

        % Continuous elapsed years from baseline visit
        matched_time(i) = days(date_current - date_base) / 365.25;
        clinical_mask(i) = true;
    else
        warning('Missing baseline or current clinical data for %s at Visit %d.', raw_sub_ids{i}, current_visit);
    end
end

X_features       = cell2mat(band_features(clinical_mask));
valid_labels     = valid_labels(clinical_mask);
valid_groups     = valid_groups(clinical_mask);
matched_age_base = matched_age_base(clinical_mask);
matched_sex      = matched_sex(clinical_mask);
matched_hand     = matched_hand(clinical_mask);
matched_time     = matched_time(clinical_mask);
matched_id       = categorical(raw_sub_ids(clinical_mask));
valid_visits     = valid_visits(clinical_mask);
final_N          = sum(clinical_mask);

X_features(isnan(X_features)) = nanmedian(X_features(:));
fprintf('Retained %d datasets with complete clinical data.\n', final_N);

n_chans = 128;
n_freqs = length(sub_freqs);
fprintf('Feature space per subject: %d channels x %d frequencies = %d total features.\n', ...
    n_chans, n_freqs, size(X_features, 2));

% -------------------------------------------------------------------------
% 4. Deterministic Continuous Spatio-Spectral NMF (via NNDSVD)
% -------------------------------------------------------------------------
X_nmf = X_features;
X_nmf(X_nmf < 0) = 0; % Enforce non-negativity boundary

k_nmf = 5; % Number of continuous spatio-spectral modules
fprintf('Computing deterministic NNDSVD initialisation for K = %d components...\n', k_nmf);

% 1. Compute deterministic starting matrices via NNDSVD
[W0, H0] = nndsvd_init(X_nmf, k_nmf);

% 2. Fit NMF using the deterministic seed (Replicates is set to 1)
fprintf('Running Spatio-Spectral NMF (ALS)...\n');
opts = statset('MaxIter', 1000, 'Display', 'final');

[W, H] = nnmf(X_nmf, k_nmf, ...
    'W0', W0, ...
    'H0', H0, ...
    'Algorithm', 'als', ...
    'Options', opts);

% 3. Resolve scaling ambiguity: unit L2-norm on H rows, absorb scale into W
norm_H = sqrt(sum(H.^2, 2));
norm_H(norm_H == 0) = eps;
H = H ./ norm_H;
W = W .* norm_H';

% 4. Deterministic sorting: sort components by mean subject expression
[~, sort_idx] = sort(mean(W, 1), 'descend');
W = W(:, sort_idx);
H = H(sort_idx, :);

% -------------------------------------------------------------------------
% 5. Normative LME Modeling with Random Slopes
% -------------------------------------------------------------------------
fprintf('Fitting Normative LME Models on %d NMF Components...\n', k_nmf);

nmf_residuals = zeros(final_N, k_nmf);
is_control = contains(valid_groups, 'CONTROL');

for k = 1:k_nmf
    lme_data = table(W(:, k), matched_age_base, matched_sex, matched_hand, matched_time, matched_id, ...
        'VariableNames', {'NMF_Score', 'AGE_baseline', 'SEX', 'HAND', 'TIME', 'ID'});

    % Train normative model exclusively on Healthy Controls using REML
    ctrl_data = lme_data(is_control, :);
    mdl = fitlme(ctrl_data, 'NMF_Score ~ AGE_baseline + SEX + TIME + (1 + TIME|ID)', 'FitMethod', 'REML');

    % Predict for ALL subjects using Fixed Effects ONLY
    expected_score = predict(mdl, lme_data, 'Conditional', false);

    % Pathological divergence from normative healthy aging
    nmf_residuals(:, k) = lme_data.NMF_Score - expected_score;
end

% -------------------------------------------------------------------------
% Visualisation 1: Spatio-Spectral Component Profiles (Topomap + Spectrum)
% -------------------------------------------------------------------------
num_comps_to_plot = min(12, k_nmf);
fh_topo = figure('Name', 'NMF Spatio-Spectral Profiles', 'Color', 'w', ...
    'Position', [100, 100, 1000, 260 * num_comps_to_plot]);

try
    cmap_seq = brewermap([], 'YlGnBu');
catch
    cmap_seq = parula(256);
end

for k = 1:num_comps_to_plot
    % Reshape component loading row back into (Channels x Frequencies)
    comp_matrix = reshape(H(k, :), [n_chans, n_freqs]);

    % Separate spatial footprint from continuous spectral tuning via rank-1 SVD
    [u, s, v] = svds(comp_matrix, 1);
    spatial_topo    = abs(u * sqrt(s(1, 1)));
    spectral_tuning = abs(v' * sqrt(s(1, 1)));

    % Subplot Left: Topography
    sh_map = subplot(num_comps_to_plot, 2, (k - 1) * 2 + 1);
    mytopoplot(spatial_topo, [], sprintf('Comp %d Topography', k), sh_map);
    max_val = max(spatial_topo);
    if max_val > 0, clim(sh_map, [0, max_val]); end
    colormap(sh_map, cmap_seq);
    colorbar(sh_map);

    % Subplot Right: Continuous Spectral Tuning Curve
    sh_spec = subplot(num_comps_to_plot, 2, (k - 1) * 2 + 2);
    plot(sh_spec, sub_freqs, spectral_tuning, 'LineWidth', 2.2, 'Color', [0.15, 0.40, 0.70]);
    grid(sh_spec, 'on'); box(sh_spec, 'on');
    xlabel(sh_spec, 'Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(sh_spec, 'Weight', 'FontSize', 10, 'FontWeight', 'bold');
    title(sh_spec, sprintf('Comp %d Spectral Tuning (%.1f-%.1f Hz)', k, target_band(1), target_band(2)), ...
        'FontSize', 11, 'FontWeight', 'bold');
    xlim(sh_spec, [target_band(1), target_band(2)]);
    set(sh_spec, 'FontSize', 9, 'LineWidth', 1.0);
end

print(fh_topo, fullfile(output_dir, sprintf('Dashboard_NMF_Profiles_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 2: Longitudinal Residual Trajectories
% -------------------------------------------------------------------------
time_months = matched_time * 12;
unique_ids  = unique(matched_id);

c_ctrl = [0.15, 0.65, 0.25];
c_als  = [0.85, 0.20, 0.20];

n_cols = min(num_comps_to_plot, 6);
n_rows = ceil(num_comps_to_plot / n_cols);

fh_traj = figure('Name', 'Longitudinal NMF Residual Trajectories', 'Color', 'w', ...
    'Position', [150, 150, 420 * n_cols, 360 * n_rows]);

max_t  = ceil(max(time_months));
t_grid = linspace(0, max_t, 100)';

for k = 1:num_comps_to_plot
    sh = subplot(n_rows, n_cols, k);
    hold(sh, 'on');
    grid(sh, 'on');

    comp_res = nmf_residuals(:, k);
    yline(sh, 0, '--k', 'LineWidth', 1.3, 'Alpha', 0.7, 'DisplayName', 'Normative Baseline');

    % Spaghetti lines per participant
    for i = 1:length(unique_ids)
        sub_mask = (matched_id == unique_ids(i));
        t_sub    = time_months(sub_mask);
        r_sub    = comp_res(sub_mask);

        [t_sub, sort_order] = sort(t_sub);
        r_sub = r_sub(sort_order);

        if is_control(find(sub_mask, 1))
            plot(sh, t_sub, r_sub, '-', 'Color', [c_ctrl, 0.20], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        else
            plot(sh, t_sub, r_sub, '-', 'Color', [c_als, 0.20], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        end
    end

    scatter(sh, time_months(is_control), comp_res(is_control), 25, c_ctrl, 'filled', ...
        'MarkerFaceAlpha', 0.55, 'DisplayName', 'Control');
    scatter(sh, time_months(~is_control), comp_res(~is_control), 25, c_als, 'filled', ...
        'MarkerFaceAlpha', 0.55, 'DisplayName', 'ALS');

    p_ctrl     = polyfit(time_months(is_control), comp_res(is_control), 1);
    y_ctrl_fit = polyval(p_ctrl, t_grid);
    plot(sh, t_grid, y_ctrl_fit, '-', 'Color', [0.05, 0.40, 0.15], 'LineWidth', 2.2, 'DisplayName', 'Control Trend');

    p_als     = polyfit(time_months(~is_control), comp_res(~is_control), 1);
    y_als_fit = polyval(p_als, t_grid);
    plot(sh, t_grid, y_als_fit, '-', 'Color', [0.55, 0.05, 0.10], 'LineWidth', 2.2, 'DisplayName', 'ALS Fit');

    title(sh, sprintf('NMF Component %d', k), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel(sh, 'Elapsed Time (Months)', 'FontSize', 10);
    ylabel(sh, 'Normative Residual', 'FontSize', 10);
    set(sh, 'FontSize', 9, 'LineWidth', 1.0, 'Box', 'on');
    xlim(sh, [-1, max_t + 2]);

    if k == 1, legend(sh, 'Location', 'best', 'FontSize', 8); end
end

print(fh_traj, fullfile(output_dir, sprintf('Dashboard_NMF_Spaghetti_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% 6. UMAP on Normative Residuals
% -------------------------------------------------------------------------
n_neighbors = min(25, floor(final_N / 3));
fprintf('Running UMAP on Normative NMF Residuals...\n');

[Y_umap, ~] = run_umap(nmf_residuals, ...
    'n_neighbors', n_neighbors, ...
    'min_dist', 0.1, ...
    'metric', 'euclidean', ...
    'method', 'MATLAB', ...
    'ask_to_save', false, ...
    'verbose', 'text');

% -------------------------------------------------------------------------
% 7. Dynamic Colour Palette & UMAP Dashboard
% -------------------------------------------------------------------------
unique_labels = unique(valid_labels);
num_categories = length(unique_labels);
colour_map = zeros(num_categories, 3);

ctrl_visits = sort(unique(valid_visits(contains(valid_groups, 'CONTROL'))));
als_visits  = sort(unique(valid_visits(~contains(valid_groups, 'CONTROL'))));

for k = 1:num_categories
    lbl = unique_labels{k};
    v_val = str2double(regexp(lbl, '\d+', 'match', 'once'));

    if contains(lbl, 'CONTROL')
        v_idx = find(ctrl_visits == v_val, 1);
        n_v = max(length(ctrl_visits), 1);
        frac = (n_v == 1) * 0.4 + (n_v > 1) * ((v_idx - 1) / (n_v - 1));
        colour_map(k, :) = (1 - frac) * [0.65, 0.90, 0.65] + frac * [0.05, 0.40, 0.15];
    else
        v_idx = find(als_visits == v_val, 1);
        n_v = max(length(als_visits), 1);
        frac = (n_v == 1) * 0.4 + (n_v > 1) * ((v_idx - 1) / (n_v - 1));
        colour_map(k, :) = (1 - frac) * [0.95, 0.60, 0.60] + frac * [0.50, 0.05, 0.10];
    end
end

fh_umap = figure('Name', 'NMF LME Residual UMAP', 'Color', 'w', 'Position', [200, 200, 900, 700]);
hold on; grid on;

for k = 1:num_categories
    idx = strcmp(valid_labels, unique_labels{k});
    scatter(Y_umap(idx, 1), Y_umap(idx, 2), 70, colour_map(k, :), 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceAlpha', 0.85, ...
        'DisplayName', unique_labels{k});
end

xlabel('UMAP Dimension 1', 'FontWeight', 'bold');
ylabel('UMAP Dimension 2', 'FontWeight', 'bold');
title(sprintf('NMF LME Residual UMAP (%.1f-%.1f Hz, N=%d)', target_band(1), target_band(2), final_N), 'FontSize', 13);
legend('Location', 'best', 'Interpreter', 'none');

plotX = 22; plotY = 16;
set(fh_umap, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh_umap, fullfile(output_dir, sprintf('Dashboard_NMF_LME_UMAP_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

fprintf('Finished. All spatio-spectral NMF dashboards saved to %s\n', output_dir);
end

function [W0, H0] = nndsvd_init(A, k)
% =========================================================================
% Non-Negative Double Singular Value Decomposition (NNDSVD) Initialisation
% Boutsidis & Gallopoulos (2008)
% Guarantees deterministic, non-negative starting factors without random seeds
% =========================================================================
[U, S, V] = svds(A, k);
W0 = zeros(size(U));
H0 = zeros(k, size(V, 1));

% First singular triplet
W0(:, 1) = sqrt(S(1, 1)) * abs(U(:, 1));
H0(1, :) = sqrt(S(1, 1)) * abs(V(:, 1))';

% Remaining k-1 components
for j = 2:k
    u = U(:, j);
    v = V(:, j);

    u_p = max(u, 0); u_n = abs(min(u, 0));
    v_p = max(v, 0); v_n = abs(min(v, 0));

    n_p = norm(u_p) * norm(v_p);
    n_n = norm(u_n) * norm(v_n);

    if n_p >= n_n
        W0(:, j) = sqrt(S(j, j) * n_p) * (u_p / max(norm(u_p), eps));
        H0(j, :) = sqrt(S(j, j) * n_p) * (v_p' / max(norm(v_p), eps));
    else
        W0(:, j) = sqrt(S(j, j) * n_n) * (u_n / max(norm(u_n), eps));
        H0(j, :) = sqrt(S(j, j) * n_n) * (v_n' / max(norm(v_n), eps));
    end
end
end