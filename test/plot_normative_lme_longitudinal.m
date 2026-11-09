function plot_normative_lme_longitudinal(base_paths, output_dir, table_path, target_band)
% =========================================================================
% Continuous Spatio-Spectral NMF & Longitudinal LME Modelling Pipeline
% 1. Extracts continuous relative power (Channels x Frequencies).
% 2. Fits deterministic NNDSVD-NMF to derive non-negative spectral modules.
% 3. Fits longitudinal LMEs testing GROUP * TIME interactions per component.
% 4. Plots individual trajectories and model-predicted marginal regression fits.
% =========================================================================

if nargin < 4 || isempty(target_band)
    target_band = [2, 30]; % Continuous physiological window
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
        f_all = freqs >= 1 & freqs <= 45;

        if isempty(sub_freqs)
            sub_freqs = freqs(f_idx);
        end

        psd_band = psd_post(:, f_idx);
        total_power = sum(psd_post(:, f_all), 2);
        total_power(total_power == 0) = eps;

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

% Cohort categorical variable with CONTROL as reference level
group_cohort = repmat({'ALS'}, final_N, 1);
group_cohort(contains(valid_groups, 'CONTROL')) = {'CONTROL'};
group_cohort = reordercats(categorical(group_cohort), {'CONTROL', 'ALS'});

X_features(isnan(X_features)) = nanmedian(X_features(:));
fprintf('Retained %d datasets with complete clinical data.\n', final_N);

n_chans = 128;
n_freqs = length(sub_freqs);
fprintf('Feature space per subject: %d channels x %d frequencies = %d total features.\n', ...
    n_chans, n_freqs, size(X_features, 2));

% =========================================================================
% Diagnostic Sweep: NMF Rank Selection (Elbow & Split-Half Stability)
% =========================================================================
k_candidates = 2:10;
n_splits     = 10; % Number of cross-validation split iterations
n_k          = length(k_candidates);

err_curve    = zeros(n_k, 1);
stability    = zeros(n_k, 1);

X_nmf = X_features;
X_nmf(X_nmf < 0) = 0;

algo = 'als'; % als / mult
opts = statset('MaxIter', 2000, 'Display', 'off');
f_norm_total = norm(X_nmf, 'fro');

fprintf('Evaluating NMF candidate ranks K = %d to %d...\n', k_candidates(1), k_candidates(end));

for ik = 1:n_k
    k = k_candidates(ik);
    fprintf('  Testing K = %d... ', k);

    % 1. Evaluate full-dataset reconstruction error via deterministic NNDSVD
    [W0, H0] = nndsvda_init(X_nmf, k);
    [W_full, H_full] = nnmf(X_nmf, k, 'W0', W0, 'H0', H0, 'Algorithm', algo, 'Options', opts);
    err_curve(ik) = norm(X_nmf - W_full * H_full, 'fro') / f_norm_total;

    % 2. Split-half stability cross-validation
    sim_scores = zeros(n_splits, 1);
    n_total    = size(X_nmf, 1);
    half_n     = floor(n_total / 2);

    for s = 1:n_splits
        perm = randperm(n_total);
        idx_A = perm(1:half_n);
        idx_B = perm(half_n+1 : 2*half_n);

        [W0_A, H0_A] = nndsvda_init(X_nmf(idx_A, :), k);
        [~, H_A]     = nnmf(X_nmf(idx_A, :), k, 'W0', W0_A, 'H0', H0_A, 'Algorithm', algo, 'Options', opts);

        [W0_B, H0_B] = nndsvda_init(X_nmf(idx_B, :), k);
        [~, H_B]     = nnmf(X_nmf(idx_B, :), k, 'W0', W0_B, 'H0', H0_B, 'Algorithm', algo, 'Options', opts);

        % Normalize component basis rows
        H_A_norm = H_A ./ sqrt(sum(H_A.^2, 2));
        H_B_norm = H_B ./ sqrt(sum(H_B.^2, 2));

        % Compute pairwise cosine similarity matrix [k x k]
        cos_sim = H_A_norm * H_B_norm';

        % Hungarian / Greedy match best-paired components
        match_sim = zeros(k, 1);
        temp_sim = cos_sim;
        for m = 1:k
            [max_val, lin_idx] = max(temp_sim(:));
            [row, col] = ind2sub([k, k], lin_idx);
            match_sim(m) = max_val;
            temp_sim(row, :) = -inf;
            temp_sim(:, col) = -inf;
        end
        sim_scores(s) = mean(match_sim);
    end

    stability(ik) = mean(sim_scores);
    fprintf('Error = %.4f | Split Stability = %.3f\n', err_curve(ik), stability(ik));
end

% -------------------------------------------------------------------------
% Plot Diagnostic Curves
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [200, 200, 900, 400]);

% Left: Reconstruction Error (Elbow Curve)
subplot(1, 2, 1);
plot(k_candidates, err_curve, '-o', 'Color', [0.15, 0.40, 0.70], ...
    'LineWidth', 2, 'MarkerFaceColor', [0.15, 0.40, 0.70]);
grid on; box on;
xlabel('Number of Components (K)', 'FontWeight', 'bold');
ylabel('Relative Reconstruction Error', 'FontWeight', 'bold');
title('Reconstruction Error (Elbow)', 'FontSize', 11);
set(gca, 'XTick', k_candidates, 'FontSize', 10);

% Right: Cross-Split Stability
subplot(1, 2, 2);
plot(k_candidates, stability, '-s', 'Color', [0.75, 0.25, 0.15], ...
    'LineWidth', 2, 'MarkerFaceColor', [0.75, 0.25, 0.15]);
yline(0.85, '--k', 'LineWidth', 1.2, 'DisplayName', 'Robust Threshold (0.85)');
grid on; box on;
xlabel('Number of Components (K)', 'FontWeight', 'bold');
ylabel('Mean Cosine Similarity', 'FontWeight', 'bold');
title('Split-Half Component Stability', 'FontSize', 11);
ylim([0.4, 1.0]);
set(gca, 'XTick', k_candidates, 'FontSize', 10);
legend('Stability across splits', 'Threshold', 'Location', 'southwest');

% -------------------------------------------------------------------------
% 4. Deterministic Continuous Spatio-Spectral NMF (via NNDSVD)
% -------------------------------------------------------------------------
X_nmf = X_features;
X_nmf(X_nmf < 0) = 0;

k_nmf = 4;
fprintf('Computing deterministic NNDSVD initialisation for K = %d components...\n', k_nmf);

[W0, H0] = nndsvda_init(X_nmf, k_nmf);

fprintf('Running Spatio-Spectral NMF (ALS)...\n');
opts = statset('MaxIter', 2000, 'Display', 'final');

[W, H] = nnmf(X_nmf, k_nmf, 'W0', W0, 'H0', H0, 'Algorithm', algo, 'Options', opts);

norm_H = sqrt(sum(H.^2, 2));
norm_H(norm_H == 0) = eps;
H = H ./ norm_H;
W = W .* norm_H';

[~, sort_idx] = sort(mean(W, 1), 'descend');
W = W(:, sort_idx);
H = H(sort_idx, :);

% -------------------------------------------------------------------------
% 5. Longitudinal LME Modelling (ALS vs CONTROL Cohort Trajectories)
% -------------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf('Fitting Longitudinal LMEs: Score ~ GROUP * TIME + AGE + SEX + HAND\n');
fprintf('==================================================================\n');

lme_models = cell(k_nmf, 1);
p_interaction = zeros(k_nmf, 1);
beta_interaction = zeros(k_nmf, 1);

for k = 1:k_nmf
    lme_data = table(W(:, k), group_cohort, matched_age_base, matched_sex, matched_hand, matched_time, matched_id, ...
        'VariableNames', {'NMF_Score', 'GROUP', 'AGE_baseline', 'SEX', 'HAND', 'TIME', 'ID'});

    % Fit combined cohort model using REML
    mdl = fitlme(lme_data, ...
        'NMF_Score ~ 1 + GROUP*TIME + AGE_baseline + SEX + (1 + TIME | ID)', ...
        'FitMethod', 'REML');

    lme_models{k} = mdl;

    % Extract GROUP_ALS:TIME interaction stats
    fe_names = mdl.CoefficientNames;
    int_idx = find(strcmp(fe_names, 'GROUP_ALS:TIME'));
    if ~isempty(int_idx)
        beta_interaction(k) = mdl.Coefficients.Estimate(int_idx);
        p_interaction(k)    = mdl.Coefficients.pValue(int_idx);
    else
        p_interaction(k) = NaN;
    end

    fprintf('Comp %2d: GROUP_ALS:TIME beta = %+7.4f, t(%d) = %+6.2f, p = %.4f\n', ...
        k, beta_interaction(k), mdl.Coefficients.DF(int_idx), ...
        mdl.Coefficients.tStat(int_idx), p_interaction(k));
end

% -------------------------------------------------------------------------
% Visualisation 1: Spatio-Spectral Component Profiles (Topomap + Spectrum)
% -------------------------------------------------------------------------
num_comps_to_plot = k_nmf;
fh_topo = figure('Name', 'NMF Spatio-Spectral Profiles', 'Color', 'w', ...
    'Position', [100, 100, 1000, 240 * num_comps_to_plot]);

try
    cmap_seq = brewermap([], 'YlGnBu');
catch
    cmap_seq = parula(256);
end

for k = 1:num_comps_to_plot
    comp_matrix = reshape(H(k, :), [n_chans, n_freqs]);

    [u, s, v] = svds(comp_matrix, 1);
    spatial_topo    = abs(u * sqrt(s(1, 1)));
    spectral_tuning = abs(v' * sqrt(s(1, 1)));

    sh_map = subplot(num_comps_to_plot, 2, (k - 1) * 2 + 1);
    mytopoplot(spatial_topo, [], sprintf('Comp %d Scalp Topography', k), sh_map);
    max_val = max(spatial_topo);
    if max_val > 0, clim(sh_map, [0, max_val]); end
    colormap(sh_map, cmap_seq);
    colorbar(sh_map);

    sh_spec = subplot(num_comps_to_plot, 2, (k - 1) * 2 + 2);
    plot(sh_spec, sub_freqs, spectral_tuning, 'LineWidth', 2.2, 'Color', [0.15, 0.40, 0.70]);
    grid(sh_spec, 'on'); box(sh_spec, 'on');
    xlabel(sh_spec, 'Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(sh_spec, 'Loading Weight', 'FontSize', 10, 'FontWeight', 'bold');
    title(sh_spec, sprintf('Comp %d Spectral Tuning (p_{int} = %.4f)', k, p_interaction(k)), ...
        'FontSize', 11, 'FontWeight', 'bold');
    xlim(sh_spec, [target_band(1), target_band(2)]);
    set(sh_spec, 'FontSize', 9, 'LineWidth', 1.0);
end

print(fh_topo, fullfile(output_dir, sprintf('Dashboard_NMF_Profiles_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 2: Longitudinal Cohort Trajectories (Data + Model Fits)
% -------------------------------------------------------------------------
time_months = matched_time * 12;
unique_ids  = unique(matched_id);
is_control  = (group_cohort == 'CONTROL');

c_ctrl = [0.15, 0.65, 0.25]; % Control: Green
c_als  = [0.85, 0.20, 0.20]; % ALS: Red

n_cols = min(num_comps_to_plot, 3);
n_rows = ceil(num_comps_to_plot / n_cols);

fh_traj = figure('Name', 'Longitudinal NMF Cohort Trajectories', 'Color', 'w', ...
    'Position', [120, 80, 430 * n_cols, 360 * n_rows]);

max_m = ceil(max(time_months));
grid_months = linspace(0, max_m, 100)';
grid_years  = grid_months / 12;

mean_age  = mean(matched_age_base);
ref_sex   = categorical({'M'}, categories(matched_sex));
ref_hand  = categorical({'R'}, categories(matched_hand));

for k = 1:num_comps_to_plot
    sh = subplot(n_rows, n_cols, k);
    hold(sh, 'on');
    grid(sh, 'on');

    comp_scores = W(:, k);
    mdl = lme_models{k};

    % 1. Individual subject trajectories (Spaghetti lines)
    for i = 1:length(unique_ids)
        sub_mask = (matched_id == unique_ids(i));
        t_sub    = time_months(sub_mask);
        s_sub    = comp_scores(sub_mask);

        [t_sub, sort_order] = sort(t_sub);
        s_sub = s_sub(sort_order);

        if is_control(find(sub_mask, 1))
            plot(sh, t_sub, s_sub, '-', 'Color', [c_ctrl, 0.20], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        else
            plot(sh, t_sub, s_sub, '-', 'Color', [c_als, 0.20], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
    end

    % 2. Individual visit observations
    scatter(sh, time_months(is_control), comp_scores(is_control), 20, c_ctrl, 'filled', ...
        'MarkerFaceAlpha', 0.50, 'DisplayName', 'Control Obs');
    scatter(sh, time_months(~is_control), comp_scores(~is_control), 20, c_als, 'filled', ...
        'MarkerFaceAlpha', 0.50, 'DisplayName', 'ALS Obs');

    % 3. Marginal population-level LME predictions (Fixed Effects only)
    % Dummy ID using the first subject category to satisfy validateDataset
    dummy_id = repmat(matched_id(1), 100, 1);

    % tbl_pred_ctrl = table(categorical(repmat({'CONTROL'}, 100, 1), categories(group_cohort)), ...
    %     repmat(mean_age, 100, 1), repmat(ref_sex, 100, 1), repmat(ref_hand, 100, 1), grid_years, dummy_id, ...
    %     'VariableNames', {'GROUP', 'AGE_baseline', 'SEX', 'HAND', 'TIME', 'ID'});
    %
    % tbl_pred_als = table(categorical(repmat({'ALS'}, 100, 1), categories(group_cohort)), ...
    %     repmat(mean_age, 100, 1), repmat(ref_sex, 100, 1), repmat(ref_hand, 100, 1), grid_years, dummy_id, ...
    %     'VariableNames', {'GROUP', 'AGE_baseline', 'SEX', 'HAND', 'TIME', 'ID'});

    tbl_pred_ctrl = table(categorical(repmat({'CONTROL'}, 100, 1), categories(group_cohort)), ...
        repmat(mean_age, 100, 1), repmat(ref_sex, 100, 1), grid_years, dummy_id, ...
        'VariableNames', {'GROUP', 'AGE_baseline', 'SEX', 'TIME', 'ID'});

    tbl_pred_als = table(categorical(repmat({'ALS'}, 100, 1), categories(group_cohort)), ...
        repmat(mean_age, 100, 1), repmat(ref_sex, 100, 1), grid_years, dummy_id, ...
        'VariableNames', {'GROUP', 'AGE_baseline', 'SEX', 'TIME', 'ID'});

    y_hat_ctrl = predict(mdl, tbl_pred_ctrl, 'Conditional', false);
    y_hat_als  = predict(mdl, tbl_pred_als,  'Conditional', false);

    plot(sh, grid_months, y_hat_ctrl, '-', 'Color', [0.05, 0.40, 0.15], 'LineWidth', 2.5, ...
        'DisplayName', 'Control Model Fit');
    plot(sh, grid_months, y_hat_als,  '-', 'Color', [0.60, 0.05, 0.10], 'LineWidth', 2.5, ...
        'DisplayName', 'ALS Model Fit');

    % Title with interaction statistical significance
    if p_interaction(k) < 0.001
        p_str = 'p < 0.001***';
    elseif p_interaction(k) < 0.05
        p_str = sprintf('p = %.3f*', p_interaction(k));
    else
        p_str = sprintf('p = %.3f (ns)', p_interaction(k));
    end

    title(sh, sprintf('Comp %d (Rate Diff: %s)', k, p_str), 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(sh, 'Elapsed Time (Months)', 'FontSize', 9);
    ylabel(sh, 'NMF Activation Score', 'FontSize', 9);
    set(sh, 'FontSize', 9, 'LineWidth', 1.0, 'Box', 'on');
    xlim(sh, [-1, max_m + 2]);

    if k == 1, legend(sh, 'Location', 'best', 'FontSize', 8); end
end

print(fh_traj, fullfile(output_dir, sprintf('Dashboard_LME_Longitudinal_Trajectories_%.1f_%.1fHz', target_band(1), target_band(2))), '-dtiff', '-r300');

fprintf('\nPipeline complete. Figures saved to: %s\n', output_dir);
end

% -------------------------------------------------------------------------
% Deterministic NNDSVD Initialisation Helper
% -------------------------------------------------------------------------
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

function [W0, H0] = nndsvda_init(A, k)
% =========================================================================
% NNDSVDa: Non-Negative Double SVD with Average Infill (No Zeros)
% Replaces zero entries with the mean of matrix A.
% Fully deterministic, strictly positive (W0 > 0, H0 > 0).
% =========================================================================
[m, n] = size(A);
[U, S, V] = svds(A, k);
W0 = zeros(m, k);
H0 = zeros(k, n);

avg_A = mean(A(:));

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

% Replace all zeros with average of A
W0(W0 == 0) = avg_A;
H0(H0 == 0) = avg_A;
end