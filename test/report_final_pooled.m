function [DATA, flag_redo] = report_final_pooled(base_paths, output_dir)
% =========================================================================
% Script for reporting on pooled EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
% =========================================================================
% INPUTS:
%   base_paths - Cell array of strings to cohort directories containing
%                subject subfolders.
%                Example: {'F:\3_PREPROCESSED_DATA\RS\ALS\T5',
%                          'F:\3_PREPROCESSED_DATA\RS\CONTROL\T1'}
%   output_dir - String specifying where to save the pooled dashboards.
% =========================================================================

fprintf('\n==================================================================\n');
fprintf('Generating the final pooled quality dashboard.\n');
fprintf('==================================================================\n');

if exist(output_dir, 'dir') ~= 7, mkdir(output_dir); end

% 1. Identify all subject folders across all provided base paths
all_subjects = {};
all_paths    = {};
cohort_names = {};

for p = 1:length(base_paths)
    current_path = regexprep(base_paths{p}, '[\\/]+$', ''); % Remove trailing slashes

    % Extract cohort from parent directory (e.g. ALS) and visit (e.g. T5)
    [parent_dir, visit_name] = fileparts(current_path);
    [~, cohort_main] = fileparts(parent_dir);

    % Combine to make 'ALS_T5'
    c_name = sprintf('%s_%s', cohort_main, visit_name);

    d = dir(current_path);
    isub = [d(:).isdir];
    subj_folders = {d(isub).name}';

    % Remove generic MATLAB directories and report folders
    subj_folders(ismember(subj_folders, {'.', '..', 'reports'})) = [];

    for s = 1:length(subj_folders)
        all_subjects{end+1} = subj_folders{s};
        all_paths{end+1}    = current_path;
        cohort_names{end+1} = c_name;
    end
end

num_subject = length(all_subjects);
if num_subject == 0
    error('No subject subfolders found in the provided base paths.');
end

% Preallocate metric arrays
valid_subjects   = false(num_subject, 1);
pooled_subjects  = cell(num_subject, 1);
cohort_labels    = cell(num_subject, 1);
N_trials_start   = NaN(num_subject, 1);
N_trials_removed = NaN(num_subject, 1);
N_interp_chan    = NaN(num_subject, 1);
N_interp_trl     = NaN(num_subject, 1);
P_extreme        = NaN(num_subject, 1);
P_emg_left       = NaN(num_subject, 1);
V_eog_left       = NaN(num_subject, 1);
V_shift          = NaN(128, num_subject);
auto_status      = NaN(num_subject, 1);
ica_removed      = NaN(num_subject, 1);
cluster_sizes    = NaN(num_subject, 1);
chan_offsets     = NaN(num_subject, 1);
cov_matrices     = NaN(128, 128, num_subject);

fprintf('Loading QA metrics for %d datasets across %d folders... \n', num_subject, length(base_paths));

for i_subject = 1:num_subject
    sub_id       = all_subjects{i_subject};
    current_path = all_paths{i_subject};
    c_name       = cohort_names{i_subject};

    % Attempt to construct a minimal myPaths for your existing pipeline function
    tmp_myPaths = struct('preproc', current_path, 'group', c_name);

    try
        subject = preproc_folders_subject(sub_id, tmp_myPaths);
        qa_name = fullfile(subject.qa, subject.qametrics{1});
    catch
        % Fallback: Dynamically search for the QA mat file if preproc_folders_subject fails
        qa_files = dir(fullfile(current_path, sub_id, '**', '*qametrics*.mat'));
        if isempty(qa_files)
            qa_files = dir(fullfile(current_path, sub_id, '**', '*qa*.mat'));
        end
        if ~isempty(qa_files)
            qa_name = fullfile(qa_files(1).folder, qa_files(1).name);
        else
            continue; % Skip if no QA file exists
        end
    end

    if exist(qa_name, 'file')
        load(qa_name, 'qa_data');
        qa_data = report_issues(qa_data);

        % Prefix the subject ID with the cohort name for plotting
        pooled_subjects{i_subject} = sprintf('%s_%s', c_name, sub_id);
        cohort_labels{i_subject}   = c_name; % Keep track for clustering plot

        % Retained and Removed Trial Counts
        N_trials_start(i_subject)   = qa_data.issues_to_check.NumberTrials1;
        N_trials_removed(i_subject) = qa_data.issues_to_check.NumberTrials1 - qa_data.issues_to_check.NumberTrials3;

        % Bad channels and cluster sizes
        bad_chans = qa_data.bad_channel.electrodes;
        N_interp_chan(i_subject) = length(bad_chans) / 128;
        if ~isempty(bad_chans)
            [~, cluster_sizes_tmp] = find_elec_clusters(bad_chans);
            cluster_sizes(i_subject) = max(cluster_sizes_tmp);
        else
            cluster_sizes(i_subject) = 0;
        end

        % Interpolated trials
        if ~isempty(qa_data.epochRejections.interpEpochs)
            total_trials  = qa_data.issues_to_check.NumberTrials1;
            interp_epochs = qa_data.epochRejections.interpEpochs;
            N_interp_trl(i_subject) = interp_epochs / total_trials;
        else
            N_interp_trl(i_subject) = 0;
        end

        % Median voltage shift
        V_shift(:, i_subject) = qa_data.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        % Automagic status
        status_str = lower(strtrim(qa_data.automagicmetrics.status));
        switch status_str
            case 'good', auto_status(i_subject) = 0;
            case 'ok',   auto_status(i_subject) = 0.5;
            case 'bad',  auto_status(i_subject) = 1;
            otherwise,   auto_status(i_subject) = NaN;
        end

        % Residual artifacts
        P_extreme(i_subject)          = qa_data.bad_segments.perc_bad;
        P_emg_left(i_subject)         = qa_data.leftovers.muscle2;
        V_eog_left(i_subject)         = qa_data.leftovers.blink_stats.peak_post_fp;
        ica_removed(i_subject)        = sum(qa_data.ica.final.removed);
        chan_offsets(i_subject)       = mean(abs(qa_data.bad_channel.offsets.offsets_mV(1:128, :)), 'all');
        cov_matrices(:, :, i_subject) = qa_data.channelcov.cov_hf(1:128, 1:128);
        valid_subjects(i_subject)     = true;
    end
end
fprintf('Done!\n');

% Filter to valid entries
plot_subjects    = pooled_subjects(valid_subjects);
plot_cohorts     = cohort_labels(valid_subjects);
audit_paths      = all_paths(valid_subjects);
audit_subjects   = all_subjects(valid_subjects);
N_trials_start   = N_trials_start(valid_subjects);
N_trials_removed = N_trials_removed(valid_subjects);
N_trials_left    = N_trials_start - N_trials_removed;
N_interp_chan    = N_interp_chan(valid_subjects);
N_interp_trl     = N_interp_trl(valid_subjects);
P_extreme        = P_extreme(valid_subjects);
P_emg_left       = P_emg_left(valid_subjects);
V_eog_left       = V_eog_left(valid_subjects);
V_shift          = V_shift(:, valid_subjects);
auto_status      = auto_status(valid_subjects);
cluster_sizes    = cluster_sizes(valid_subjects);
chan_offsets     = chan_offsets(valid_subjects);
ica_removed      = ica_removed(valid_subjects);
cov_matrices     = cov_matrices(:, :, valid_subjects);
NSUB_valid       = length(plot_subjects);
V_shift_med      = median(V_shift, 1)';

valid_trl = ~all(isnan(N_interp_trl)) && any(N_interp_trl > 0);

% -------------------------------------------------------------------------
% Visualisation 5: PCA Clustering (2D space by Cohort)
% -------------------------------------------------------------------------
try
    % Compile all numerical QA metrics into a feature matrix
    if valid_trl
        X_feats = [N_trials_left, N_trials_removed, N_interp_chan, N_interp_trl, P_extreme, P_emg_left, V_eog_left, chan_offsets, cluster_sizes, ica_removed, V_shift_med];
    else
        X_feats = [N_trials_left, N_trials_removed, N_interp_chan, P_extreme, P_emg_left, V_eog_left, chan_offsets, cluster_sizes, ica_removed, V_shift_med];
    end

    % Replace any stray NaNs with median
    X_feats(isnan(X_feats)) = nanmedian(X_feats(:));

    % Z-score normalize features so variables with large scales (e.g. trials) don't dominate
    X_norm = (X_feats - mean(X_feats, 1)) ./ std(X_feats, 0, 1);

    % Perform PCA
    [~, score, ~, ~, explained] = pca(X_norm);

    fh5 = figure('Name', 'QA Metrics 2D Clustering', 'Color', 'w', 'Position', [200, 200, 800, 600]);
    hold on; grid on;

    unique_cohorts = unique(plot_cohorts);
    colors = lines(length(unique_cohorts));

    for c = 1:length(unique_cohorts)
        idx = strcmp(plot_cohorts, unique_cohorts{c});
        scatter(score(idx, 1), score(idx, 2), 70, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', unique_cohorts{c});
    end

    xlabel(sprintf('Principal Component 1 (%.1f%% variance)', explained(1)), 'FontWeight', 'bold');
    ylabel(sprintf('Principal Component 2 (%.1f%% variance)', explained(2)), 'FontWeight', 'bold');
    title('2D Clustering of Processing QA Metrics by Cohort/Visit', 'FontSize', 14);
    legend('Location', 'best', 'Interpreter', 'none');

    plotX = 22; plotY = 16;
    set(fh5, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
        'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
    print(fh5, fullfile(output_dir, 'Dashboard_PCA_Clustering'), '-dtiff', '-r300');
catch ME
    fprintf('Could not generate PCA clustering plot: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Visualisation 4: Channel Covariance Matrix Check
% -------------------------------------------------------------------------
[deviant_cov_indices, fh4] = check_channelcov(cov_matrices, plot_subjects);
plotX = 35; plotY = 45;
set(fh4, 'InvertHardCopy', 'Off', 'Color', [1 1 1]);
set(fh4, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
    'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh4, fullfile(output_dir, 'Summary_Pooled_Covariance'), '-dtiff', '-r400');

% -------------------------------------------------------------------------
% Outlier Determination for Figure 1
% -------------------------------------------------------------------------
flag_trials_lost = flag_upper_outlier(N_trials_removed, 1);
flag_trials_low  = flag_lower_outlier(N_trials_left, 30);
flag_trials      = flag_trials_lost | flag_trials_low;
flag_interp_chan = flag_upper_outlier(N_interp_chan, 0.10);
flag_cluster     = flag_upper_outlier(cluster_sizes, 3);
flag_extreme     = flag_upper_outlier(P_extreme, 0.05);
flag_emg         = flag_upper_outlier(P_emg_left);
flag_eog         = flag_upper_outlier(V_eog_left);
flag_offset      = flag_upper_outlier(chan_offsets);
flag_cov         = false(NSUB_valid, 1);
if ~isempty(deviant_cov_indices)
    flag_cov(deviant_cov_indices) = true;
end
flag_auto        = (auto_status == 1);

Flag_Matrix = [flag_trials, flag_interp_chan, flag_cluster, ...
    flag_extreme, flag_emg, flag_eog, flag_offset, flag_cov, flag_auto];
heatmap_labels = {'Trials Lost', 'Interp Chans', 'Elec Cluster', ...
    'Extreme Data', 'Leftover EMG', 'Leftover EOG', ...
    'Mean Offset', 'Covariance Dev', 'Automagic Bad'};

if valid_trl
    flag_interp_trl = flag_upper_outlier(N_interp_trl, 0.05);
    Flag_Matrix = [Flag_Matrix(:, 1:2), flag_interp_trl, Flag_Matrix(:, 3:end)];
    heatmap_labels = [heatmap_labels(1:2), {'Interp Trials'}, heatmap_labels(3:end)];
end

Total_Flags = sum(Flag_Matrix, 2);
[Total_Flags_Sorted, sort_idx] = sort(Total_Flags, 'descend');
subjects_sorted    = plot_subjects(sort_idx);
Flag_Matrix_Sorted = Flag_Matrix(sort_idx, :);

% -------------------------------------------------------------------------
% Visualisation 1: Objective Defect Co-occurrence Matrix
% -------------------------------------------------------------------------
fh1 = figure('Name', 'Quality Assurance Defect Matrix (Pooled)', 'Color', 'w', 'Position', [100, 100, 1200, 850]);
t1 = tiledlayout(1, 4, 'TileSpacing', 'compact');

nexttile([1 3]);
h = heatmap(heatmap_labels, subjects_sorted, double(Flag_Matrix_Sorted));
h.Title           = 'Objective Outlier Matrix (Pooled Cohorts)';
h.Colormap        = [0.93 0.93 0.95; 0.85 0.15 0.15];
h.ColorLimits     = [0 1];
h.ColorbarVisible = 'off';
h.XLabel          = 'Evaluated Quality Domains';
h.YLabel          = 'Participants (Most Deviant to Cleanest)';

nexttile;
b = barh(Total_Flags_Sorted, 'FaceColor', 'flat', 'EdgeColor', 'none');
for i_b = 1:length(Total_Flags_Sorted)
    if Total_Flags_Sorted(i_b) >= 2
        b.CData(i_b, :) = [0.85 0.15 0.15];
    elseif Total_Flags_Sorted(i_b) == 1
        b.CData(i_b, :) = [0.95 0.65 0.15];
    else
        b.CData(i_b, :) = [0.25 0.65 0.35];
    end
end
set(gca, 'YDir', 'reverse');
ylim([0.5, NSUB_valid + 0.5]);
yticks([]);
title('Total Red Flags');
xlabel('Defect Count');
xline(1.5, 'r--', 'Suspect (\geq 2)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

plotX = 30; plotY = max(15, NSUB_valid * 0.5);
set(fh1, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
    'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh1, fullfile(output_dir, 'Dashboard_DefectMatrix_Pooled'), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 2: Distribution Swarmcharts (4 Columns)
% -------------------------------------------------------------------------
metrics = [{N_trials_left, N_trials_removed}, ...
    {N_interp_chan, P_extreme, P_emg_left, V_eog_left, V_shift_med, cluster_sizes, ica_removed, chan_offsets}];
titles  = [{'Trials Left', 'Trials Removed'}, ...
    {'Interpolated Channels (%)', 'Extreme Data (%)', 'Leftover EMG (%)', 'Leftover EOG Peak (uV)', ...
    'Median Voltage Shift (uV)', 'Max Bad Elec Cluster Size', 'Number of ICs Removed', 'Mean Channel Offset (mV)'}];

if valid_trl
    metrics = [metrics(1:3), {N_interp_trl}, metrics(4:end)];
    titles  = [titles(1:3), {'Interpolated Trials (%)'}, titles(4:end)];
end

num_metrics = length(metrics);
Ncol = 4;
Nrow = ceil(num_metrics / Ncol);
fh2 = figure('Name', 'Pooled Cohort Distributions', 'Color', 'w', 'Position', [150, 150, Ncol*320, Nrow*260]);
t2 = tiledlayout(Nrow, Ncol, 'TileSpacing', 'loose', 'Padding', 'compact');

for i_metric = 1:num_metrics
    nexttile; hold on;
    raw_data = metrics{i_metric};
    valid_mask = ~isnan(raw_data);
    data = raw_data(valid_mask);
    subs_valid = plot_subjects(valid_mask);
    if isempty(data), continue; end

    is_trials_left = strcmp(titles{i_metric}, 'Trials Left') || strcmp(titles{i_metric}, 'Number of Trials (Left)');
    is_count_var   = strcmp(titles{i_metric}, 'Max Bad Elec Cluster Size') || ...
        strcmp(titles{i_metric}, 'Max {Bad Channels Cluster}') || ...
        strcmp(titles{i_metric}, 'Number of ICs Removed') || ...
        strcmp(titles{i_metric}, 'ICs Removed') || ...
        contains(titles{i_metric}, 'Trials');

    if ~strcmp(titles{i_metric}, 'Median Voltage Shift (uV)') && ~is_count_var
        low_s  = -eps;
        high_s = max(data) * 1.2 + eps;
        extraArgs = {'BoundaryCorrection', 'reflection', 'Support', [low_s, high_s]};
    else
        low_s  = min(data) - (std(data) * 0.5);
        high_s = max(data) + (std(data) * 0.5);
        extraArgs = {'Support', [low_s, high_s]};
    end

    if std(data) > 1e-6
        [f, x_eval] = ksdensity(data, extraArgs{:});
        f_scaled = (f / max(f)) * 0.5;
        fill(x_eval, f_scaled + 1.3, [0.2 0.4 0.7], 'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.4 0.7]);
    end

    y_base = 1.0;
    y_jitter = y_base + (rand(length(data), 1) - 0.5) * 0.3;
    scatter(data, y_jitter, 25, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.3);

    if is_trials_left
        metric_outliers = flag_lower_outlier(data, 30);
    else
        metric_outliers = flag_upper_outlier(data);
    end

    outlier_idx = find(metric_outliers);
    if ~isempty(outlier_idx)
        scatter(data(outlier_idx), y_jitter(outlier_idx), 50, 'r', 'filled', 'MarkerEdgeColor', 'k');
        for w = 1:length(outlier_idx)
            idx_w = outlier_idx(w);
            text(data(idx_w), y_jitter(idx_w) + 0.18, subs_valid{idx_w}, ...
                'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Interpreter', 'none');
        end
    end

    q = prctile(data, [25 50 75]);
    line([q(2) q(2)], [0.7 1.3], 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5);
    line([q(1) q(3)], [y_base y_base], 'Color', 'k', 'LineWidth', 1.5);
    title(titles{i_metric}, 'FontSize', 11, 'FontWeight', 'bold');

    x_pad = max(std(data) * 0.2, 1e-3);
    xlim([min(data) - x_pad, max(data) + x_pad]);
    set(gca, 'YTick', [], 'YColor', 'none', 'Box', 'off', 'TickDir', 'out');
    grid on; ax = gca; ax.GridAlpha = 0.1;
end

plotX = Ncol * 10; plotY = max(12, Nrow * 8);
set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
    'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(output_dir, 'Dashboard_Distributions_Pooled'), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 3: Dedicated Voltage Shift Audit
% -------------------------------------------------------------------------
[~, worst_volt_idx] = maxk(V_shift_med, min(3, NSUB_valid));

fprintf('\nTriggering visual audit for highest voltage shift subjects...\n');
audit_voltage_offenders_pooled(audit_subjects(worst_volt_idx), audit_paths(worst_volt_idx), plot_subjects(worst_volt_idx), output_dir);

% Outputs
DATA.subjects         = plot_subjects;
DATA.cohorts          = plot_cohorts;
DATA.flag_matrix      = Flag_Matrix;
DATA.total_flags      = Total_Flags;
DATA.flag_labels      = heatmap_labels;
DATA.suspect_subjects = plot_subjects(Total_Flags >= 2);
flag_redo             = DATA.suspect_subjects;

fprintf('Pooled dashboards and voltage audit saved successfully to %s\n', output_dir);

end

% -------------------------------------------------------------------------
% Helper: Dedicated Voltage Shift Audit (Pooled version)
% -------------------------------------------------------------------------
function audit_voltage_offenders_pooled(subject_ids, subject_paths, plot_names, output_dir)

NSUB = length(subject_ids);
if NSUB == 0, return; end

screen_size = get(0, 'ScreenSize');
monitor_h = screen_size(4);
requested_h = 320 * NSUB;
final_h = min(requested_h, monitor_h * 0.85);

fh3 = figure('Name', 'Audit: Voltage Shift Inspection (Pooled)', ...
    'Color', 'w', 'Position', [100, 80, 1350, final_h]);
t3 = tiledlayout(NSUB, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:NSUB
    sub_id    = subject_ids{i};
    base_path = subject_paths{i};
    plot_id   = plot_names{i};

    % Dynamically find the cleandata file to avoid needing myPaths.visit/task
    clean_files = dir(fullfile(base_path, sub_id, '**', '*_cleandata_b.mat'));

    if ~isempty(clean_files)
        path_data_full = fullfile(clean_files(1).folder, clean_files(1).name);
        temp = load(path_data_full, 'EEG');
        EEG = temp.EEG;

        volt_shift = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);
        [~, worst_idx] = maxk(volt_shift, 4);
        worst_labels = {EEG.chanlocs(worst_idx).labels};
        worst_data = EEG.data(worst_idx, :);

        % Subplot 1: Topoplot
        ax1 = nexttile;
        myCmap = brewermap(128, 'Reds');
        topoplot(volt_shift, EEG.chanlocs(1:128), 'maplimits', [prctile(volt_shift, 5), max(volt_shift)], ...
            'headrad', 0.5, 'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'off', ...
            'style', 'map', 'shading', 'interp', ...
            'emarker2', {worst_idx, 'o', 'k', 6, 1});
        title(sprintf('%s: Voltage Topography', plot_id), 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        hcb = colorbar; hcb.Title.String = '\muV';

        % Subplot 2: PSD
        ax2 = nexttile; hold on;
        window = EEG.srate * 2;
        noverlap = 0;
        [pxx, f] = pwelch(worst_data', window, noverlap, window, EEG.srate);
        plot(f, 10*log10(mean(pxx, 2)), 'k', 'LineWidth', 2);
        patch([8 13 13 8], [-30 -30 30 30], [0 0.4 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        title(['PSD Worst: ' strjoin(worst_labels, ', ')], 'FontSize', 9, 'FontWeight', 'bold');
        xlabel('Frequency (Hz)');
        ylabel('dB');
        xlim([1 45]);
        ylim([-30 30]);
        grid on;

        % Subplot 3: 5s Trace
        ax3 = nexttile; hold on;
        mid_point = floor(size(worst_data, 2) / 2);
        win_samples = 5 * EEG.srate;
        start_idx = max(1, mid_point - floor(win_samples/2));
        end_idx = min(size(worst_data, 2), start_idx + win_samples - 1);
        t_vec = (0:(end_idx - start_idx)) / EEG.srate;
        plot(t_vec, worst_data(:, start_idx:end_idx)', 'LineWidth', 1.1);
        title('5s Trace of Worst Channels', 'FontSize', 9, 'FontWeight', 'bold');
        xlabel('Time (s)');
        ylabel('\muV');
        ylim([-100 100]);
        grid on;
    else
        nexttile([1 3]);
        text(0.5, 0.5, sprintf('Cleandata missing for %s', plot_id), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'none');
        axis off;
    end
end

plotX = 35; plotY = max(10, NSUB * 8);
set(fh3, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh3, fullfile(output_dir, 'Dashboard_Audit_Voltage_Pooled'), '-dtiff', '-r300');

end

% -------------------------------------------------------------------------
% Helper: Upper Outlier Detection
% -------------------------------------------------------------------------
function tf = flag_upper_outlier(vec, abs_min)
if nargin < 2, abs_min = -Inf; end
tf = false(size(vec));
v_clean = vec(~isnan(vec));
if isempty(v_clean), return; end

q = prctile(v_clean, [25 75]);
iqr_val = q(2) - q(1);

if iqr_val > 1e-6
    thresh = q(2) + 1.5 * iqr_val;
else
    med_val = median(v_clean);
    mad_val = mad(v_clean, 1);
    if mad_val > 1e-6
        thresh = med_val + 3 * 1.4826 * mad_val;
    else
        thresh = med_val + 2 * std(v_clean);
    end
end

thresh = max(thresh, abs_min);
tf = (vec > thresh) & ~isnan(vec);
end

% -------------------------------------------------------------------------
% Helper: Lower Outlier Detection
% -------------------------------------------------------------------------
function tf = flag_lower_outlier(vec, abs_max)
if nargin < 2, abs_max = Inf; end
tf = false(size(vec));
v_clean = vec(~isnan(vec));
if isempty(v_clean), return; end

q = prctile(v_clean, [25 75]);
iqr_val = q(2) - q(1);

if iqr_val > 1e-6
    thresh = q(1) - 1.5 * iqr_val;
else
    med_val = median(v_clean);
    mad_val = mad(v_clean, 1);
    if mad_val > 1e-6
        thresh = med_val - 3 * 1.4826 * mad_val;
    else
        thresh = med_val - 2 * std(v_clean);
    end
end

thresh = min(thresh, abs_max);
tf = (vec < thresh) & ~isnan(vec);
end