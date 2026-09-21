function [DATA, flag_redo] = report_final(myPaths, subjects)
% =========================================================================
% Script for reporting on EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
% =========================================================================

fprintf('\n==================================================================\n');
fprintf('%s: Generating the final quality dashboard.\n', myPaths.group);
fprintf('==================================================================\n');

reports_dir = fullfile(myPaths.preproc, 'reports');
if exist(reports_dir, 'dir') ~= 7, mkdir(reports_dir); end

% Preallocate metric arrays
NSUB = length(subjects);
valid_subjects   = false(NSUB, 1);
N_trials_start   = NaN(NSUB, 1);
N_trials_removed = NaN(NSUB, 1);
N_interp_chan    = NaN(NSUB, 1);
N_interp_trl     = NaN(NSUB, 1);
P_emg_left       = NaN(NSUB, 1);
V_eog_left       = NaN(NSUB, 1);
V_shift          = NaN(128, NSUB);
gamma_cv         = NaN(NSUB, 1);
auto_status      = NaN(NSUB, 1);
ica_removed      = NaN(NSUB, 1);
cluster_sizes    = NaN(NSUB, 1);
chan_offsets     = NaN(NSUB, 1);
cov_matrices     = NaN(128, 128, NSUB);
power_diff       = NaN;

fprintf('Loading QA metrics for %d datasets... ', NSUB);
for i_subj = 1:NSUB
    i_file = 1;
    subject  = preproc_folders_subject(subjects{i_subj}, myPaths);
    qa_name  = fullfile(subject.qa, subject.qametrics{i_file});

    if exist(qa_name, 'file')
        load(qa_name, 'qa_data');
        qa_data = report_issues(qa_data);

        % Retained and Removed Trial Counts
        N_trials_start(i_subj)   = qa_data.issues_to_check.NumberTrials1;
        N_trials_removed(i_subj) = qa_data.issues_to_check.NumberTrials1 - qa_data.issues_to_check.NumberTrials3;

        % Bad channels and cluster sizes
        bad_chans = qa_data.badchaninfo.badElectrodes;
        N_interp_chan(i_subj) = length(bad_chans) / 128;
        if ~isempty(bad_chans)
            [~, cluster_sizes_tmp] = find_elec_clusters(bad_chans);
            cluster_sizes(i_subj) = max(cluster_sizes_tmp);
        else
            cluster_sizes(i_subj) = 0;
        end

        % Interpolated trials
        if ~isempty(qa_data.epochRejections.interpEpochs)
            total_trials  = qa_data.issues_to_check.NumberTrials1;
            interp_epochs = qa_data.epochRejections.interpEpochs;
            N_interp_trl(i_subj) = interp_epochs / total_trials;
        else
            N_interp_trl(i_subj) = 0;
        end

        % Median voltage shift
        V_shift(:, i_subj) = qa_data.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        % Gamma spread
        if isfield(qa_data, 'psd_gamma_cv')
            gamma_cv(i_subj) = qa_data.psd_gamma_cv;
        end

        % Automagic status
        status_str = lower(strtrim(qa_data.automagicmetrics.status));
        switch status_str
            case 'good', auto_status(i_subj) = 0;
            case 'ok',   auto_status(i_subj) = 0.5;
            case 'bad',  auto_status(i_subj) = 1;
            otherwise,   auto_status(i_subj) = NaN;
        end

        % Residual artifacts
        P_emg_left(i_subj)         = qa_data.leftovers.muscle2;
        V_eog_left(i_subj)         = qa_data.leftovers.blink_stats.peak_post_fp;
        ica_removed(i_subj)        = sum(qa_data.ica.final.removed);
        chan_offsets(i_subj)       = mean(abs(qa_data.badchaninfo.offsets.offsets_mV(1:128, :)), 'all');
        cov_matrices(:, :, i_subj) = qa_data.channelcov.cov_hf(1:128, 1:128);
        power_diff(:, :, i_subj)   = qa_data.power_diff;
        valid_subjects(i_subj)     = true;
    end
end
fprintf('Done!\n');

% Filter to valid entries
subjects         = subjects(valid_subjects);
N_trials_start   = N_trials_start(valid_subjects);
N_trials_removed = N_trials_removed(valid_subjects);
N_trials_left    = N_trials_start - N_trials_removed;
N_interp_chan    = N_interp_chan(valid_subjects);
N_interp_trl     = N_interp_trl(valid_subjects);
P_emg_left       = P_emg_left(valid_subjects);
V_eog_left       = V_eog_left(valid_subjects);
V_shift          = V_shift(:, valid_subjects);
auto_status      = auto_status(valid_subjects);
cluster_sizes    = cluster_sizes(valid_subjects);
chan_offsets     = chan_offsets(valid_subjects);
ica_removed      = ica_removed(valid_subjects);
cov_matrices     = cov_matrices(:, :, valid_subjects);
NSUB_valid       = length(subjects);
V_shift_med      = median(V_shift, 1)';

valid_trl = ~all(isnan(N_interp_trl)) && any(N_interp_trl > 0);

% -------------------------------------------------------------------------
% Visualisation 4: Channel Covariance Matrix Check
% -------------------------------------------------------------------------
[deviant_cov_indices, fh4] = check_channelcov(cov_matrices, subjects);
plotX = 35; plotY = 45;
set(fh4, 'InvertHardCopy', 'Off', 'Color', [1 1 1]);
set(fh4, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
    'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh4, fullfile(reports_dir, ['Summary3_' myPaths.group '_T' num2str(myPaths.visit) '_' myPaths.task '_' myPaths.proctime]), '-dtiff', '-r400');

% -------------------------------------------------------------------------
% Outlier Determination for Figure 1
% (Excludes raw V_shift_med to prevent penalising strong alpha generators)
% -------------------------------------------------------------------------
flag_trials_lost = flag_upper_outlier(N_trials_removed, 1);
flag_trials_low  = flag_lower_outlier(N_trials_left, 30);
flag_trials      = flag_trials_lost | flag_trials_low;
flag_interp_chan = flag_upper_outlier(N_interp_chan, 0.10);
flag_cluster     = flag_upper_outlier(cluster_sizes, 3);
flag_emg         = flag_upper_outlier(P_emg_left);
flag_eog         = flag_upper_outlier(V_eog_left);
flag_offset      = flag_upper_outlier(chan_offsets);
flag_cov         = false(NSUB_valid, 1);
if ~isempty(deviant_cov_indices)
    flag_cov(deviant_cov_indices) = true;
end
flag_auto        = (auto_status == 1);

Flag_Matrix = [flag_trials, flag_interp_chan, flag_cluster, ...
    flag_emg, flag_eog, flag_offset, flag_cov, flag_auto];
heatmap_labels = {'Trials Lost', 'Interp Chans', 'Elec Cluster', ...
    'Leftover EMG', 'Leftover EOG', ...
    'Mean Offset', 'Covariance Dev', 'Automagic Bad'};

if valid_trl
    flag_interp_trl = flag_upper_outlier(N_interp_trl, 0.05);
    Flag_Matrix = [Flag_Matrix(:, 1:2), flag_interp_trl, Flag_Matrix(:, 3:end)];
    heatmap_labels = [heatmap_labels(1:2), {'Interp Trials'}, heatmap_labels(3:end)];
end

Total_Flags = sum(Flag_Matrix, 2);
[Total_Flags_Sorted, sort_idx] = sort(Total_Flags, 'descend');
subjects_sorted    = subjects(sort_idx);
Flag_Matrix_Sorted = Flag_Matrix(sort_idx, :);

% -------------------------------------------------------------------------
% Visualisation 1: Objective Defect Co-occurrence Matrix
% -------------------------------------------------------------------------
fh1 = figure('Name', 'Quality Assurance Defect Matrix', 'Color', 'w', 'Position', [100, 100, 1200, 850]);
t1 = tiledlayout(1, 4, 'TileSpacing', 'compact');

nexttile([1 3]);
h = heatmap(heatmap_labels, subjects_sorted, double(Flag_Matrix_Sorted));
h.Title       = 'Objective Outlier Matrix (Cohort Distributions & Covariance)';
h.Colormap    = [0.93 0.93 0.95; 0.85 0.15 0.15]; % Grey = Clear, Red = Outlier
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.XLabel      = 'Evaluated Quality Domains';
h.YLabel      = 'Participants (Most Deviant to Cleanest)';

nexttile;
b = barh(Total_Flags_Sorted, 'FaceColor', 'flat', 'EdgeColor', 'none');
for i_b = 1:length(Total_Flags_Sorted)
    if Total_Flags_Sorted(i_b) >= 2
        b.CData(i_b, :) = [0.85 0.15 0.15]; % Multi-defect suspect
    elseif Total_Flags_Sorted(i_b) == 1
        b.CData(i_b, :) = [0.95 0.65 0.15]; % Single warning
    else
        b.CData(i_b, :) = [0.25 0.65 0.35]; % Clean
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
print(fh1, fullfile(reports_dir, ['Dashboard_DefectMatrix_' myPaths.group]), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 2: Distribution Swarmcharts
% -------------------------------------------------------------------------
metrics = [{N_trials_left, N_trials_removed}, ...
    {N_interp_chan, P_emg_left, V_eog_left, V_shift_med, cluster_sizes, ica_removed, chan_offsets}];
titles  = [{'Number of Trials (Left)', 'Number of Trials (Removed)'}, ...
    {'Interpolated Channels (%)', 'Leftover EMG (%)', 'Leftover EOG Peak (uV)', 'Median Voltage Shift (uV)', ...
    'Max Bad Elec Cluster Size', 'Number of ICs Removed', 'Mean Channel Offset (mV)'}];

if valid_trl
    metrics = [metrics(1:3), {N_interp_trl}, metrics(4:end)];
    titles  = [titles(1:3), {'Interpolated Trials (%)'}, titles(4:end)];
end

num_metrics = length(metrics);
Ncol = 3;
Nrow = ceil(num_metrics / Ncol);
fh2 = figure('Name', 'Cohort Distributions', 'Color', 'w', 'Position', [150, 150, Ncol*400, Nrow*320]);
t2 = tiledlayout(Nrow, Ncol, 'TileSpacing', 'loose', 'Padding', 'compact');

for i_metric = 1:num_metrics
    nexttile; hold on;
    raw_data = metrics{i_metric};
    valid_mask = ~isnan(raw_data);
    data = raw_data(valid_mask);
    subs_valid = subjects(valid_mask);
    if isempty(data), continue; end

    is_trials_left = strcmp(titles{i_metric}, 'Number of Trials (Left)');
    is_count_var   = strcmp(titles{i_metric}, 'Max Bad Elec Cluster Size') || ...
        strcmp(titles{i_metric}, 'Number of ICs Removed') || ...
        contains(titles{i_metric}, 'Number of Trials');

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
    title(titles{i_metric}, 'FontSize', 12, 'FontWeight', 'bold');

    x_pad = max(std(data) * 0.2, 1e-3);
    xlim([min(data) - x_pad, max(data) + x_pad]);
    set(gca, 'YTick', [], 'YColor', 'none', 'Box', 'off', 'TickDir', 'out');
    grid on; ax = gca; ax.GridAlpha = 0.1;
end

plotX = Ncol * 17; plotY = max(15, Nrow * 10);
set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', ...
    'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(reports_dir, ['Dashboard_Distributions_' myPaths.group]), '-dtiff', '-r300');

% -------------------------------------------------------------------------
% Visualisation 3: Dedicated Voltage Shift Audit (Alpha vs Artifact Check)
% -------------------------------------------------------------------------
[~, worst_volt_idx] = maxk(V_shift_med, min(3, NSUB_valid));
worst_volt_subjects = subjects(worst_volt_idx);
fprintf('\nTriggering visual audit for highest voltage shift subjects (Alpha vs Artifact check)...\n');
audit_voltage_offenders(myPaths, worst_volt_subjects);

% Outputs
DATA.subjects         = subjects;
DATA.flag_matrix      = Flag_Matrix;
DATA.total_flags      = Total_Flags;
DATA.flag_labels      = heatmap_labels;
DATA.suspect_subjects = subjects(Total_Flags >= 2);
flag_redo             = DATA.suspect_subjects;

fprintf('Dashboards and voltage audit saved successfully to %s\n', reports_dir);

end

% -------------------------------------------------------------------------
% Helper: Dedicated Voltage Shift Audit (Alpha vs Noise Discrimination)
% -------------------------------------------------------------------------
function audit_voltage_offenders(myPaths, flagged_subjects)

NSUB = length(flagged_subjects);
if NSUB == 0, return; end

reports_dir = fullfile(myPaths.preproc, 'reports');
screen_size = get(0, 'ScreenSize');
monitor_h = screen_size(4);
requested_h = 320 * NSUB;
final_h = min(requested_h, monitor_h * 0.85);

fh3 = figure('Name', 'Audit: Voltage Shift Inspection (Alpha vs Artifact)', ...
    'Color', 'w', 'Position', [100, 80, 1350, final_h]);
t3 = tiledlayout(NSUB, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:NSUB
    subject_id = flagged_subjects{i};
    path_data = fullfile(myPaths.preproc, subject_id, myPaths.codever, 'data');
    path_data_full = fullfile(path_data, [subject_id '_T' num2str(myPaths.visit) '_' myPaths.task '_cleandata_b.mat']);

    if exist(path_data_full, 'file') == 2
        temp = load(path_data_full, 'EEG');
        EEG = temp.EEG;

        volt_shift = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);
        [~, worst_idx] = maxk(volt_shift, 4);
        worst_labels = {EEG.chanlocs(worst_idx).labels};
        worst_data = EEG.data(worst_idx, :);

        % Subplot 1: Topoplot (Alpha = Occipital focus, Artifact = Frontal/Perimeter)
        ax1 = nexttile;
        myCmap = brewermap(128, 'Reds');
        topoplot(volt_shift, EEG.chanlocs(1:128), 'maplimits', [prctile(volt_shift, 5), max(volt_shift)], ...
            'headrad', 0.5, 'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'off', ...
            'style', 'map', 'shading', 'interp', ...
            'emarker2', {worst_idx, 'o', 'k', 6, 1});
        title(sprintf('%s: Voltage Topography', subject_id), 'FontSize', 10, 'FontWeight', 'bold');
        hcb = colorbar; hcb.Title.String = '\muV';

        % Subplot 2: PSD (Alpha = Narrow 8-12 Hz peak, Artifact = Elevated/Broadband)
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

        % Subplot 3: 5s Trace (Alpha = Sinusoidal, Artifact = Railing/Drift/Jumps)
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
        text(0.5, 0.5, sprintf('Data missing for %s', subject_id), 'HorizontalAlignment', 'center', 'FontSize', 12);
        axis off;
    end
end

plotX = 35; plotY = max(10, NSUB * 8);
set(fh3, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh3, fullfile(reports_dir, ['Dashboard_Audit_Voltage_' myPaths.group]), '-dtiff', '-r300');

end

% -------------------------------------------------------------------------
% Helper: Upper Outlier Detection (Tukey Fence / MAD Fallback)
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
% Helper: Lower Outlier Detection (Tukey Fence / MAD Fallback)
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