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
ica_redo         = NaN(NSUB, 1);
ica_removed      = NaN(NSUB, 1);
cluster_sizes    = NaN(NSUB, 1);
chan_offsets     = NaN(NSUB, 1);
cov_matrices     = NaN(128, 128, NSUB);

fprintf('Loading QA metrics for %d datasets... ', NSUB);
for i_subj = 1:NSUB
    i_file = 1;
    subject  = preproc_folders_subject(subjects{i_subj}, myPaths, 2);
    qa_name  = fullfile(subject.qa, subject.qametrics{i_file});

    if exist(qa_name, 'file')
        % Load QA struct
        load(qa_name, 'qa_data');

        % Record warnings about potential issues
        qa_data = report_issues(qa_data);

        % 0. Retained and Removed Trial Counts
        if isfield(qa_data, 'issues_to_check')
            if isfield(qa_data.issues_to_check, 'NumberTrials2')
                N_trials_start(i_subj) = qa_data.issues_to_check.NumberTrials2;
            end
            if isfield(qa_data.issues_to_check, 'NumberTrials2') && isfield(qa_data.issues_to_check, 'NumberTrials3')
                N_trials_removed(i_subj) = qa_data.issues_to_check.NumberTrials2 - qa_data.issues_to_check.NumberTrials3;
            end
        end

        % 1. Interpolated Channels (%)
        bad_chans = qa_data.badchaninfo.badElectrodes;
        N_interp_chan(i_subj) = length(bad_chans) / 128;

        % Check for electrode cluster size
        if ~isempty(bad_chans)
            [~, cluster_sizes_tmp] = find_elec_clusters(bad_chans);
            cluster_sizes(i_subj) = max(cluster_sizes_tmp);
        else
            cluster_sizes(i_subj) = 0;
        end

        % 2. Interpolated Trials (%)
        if ~isempty(qa_data.epochRejections.interpEpochs)
            total_trials  = qa_data.issues_to_check.NumberTrials1;
            interp_epochs = qa_data.epochRejections.interpEpochs;
            N_interp_trl(i_subj) = interp_epochs / total_trials;
        else
            N_interp_trl(i_subj) = 0;
        end

        % 4. Voltage Shift
        V_shift(:, i_subj) = qa_data.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        % 5. Gamma Spread
        if isfield(qa_data, 'psd_gamma_cv')
            gamma_cv(i_subj) = qa_data.psd_gamma_cv;
        end

        % 6. Automagic Status
        status_str = lower(strtrim(qa_data.automagicmetrics.status));
        switch status_str
            case 'good', auto_status(i_subj) = 0;   % Green
            case 'ok',   auto_status(i_subj) = 0.5; % Yellow
            case 'bad',  auto_status(i_subj) = 1;   % Red
            otherwise,   auto_status(i_subj) = NaN;
        end

        % 3. EMG Leftover
        P_emg_left(i_subj) = qa_data.leftovers.muscle2;

        % 7. EOG Leftover
        V_eog_left(i_subj) = qa_data.leftovers.blink2.blink_stats.peak_post_fp;
        ica_redo(i_subj) = double(qa_data.leftovers.blink1.flag_redo);

        % 8. ICs Removed
        ica_removed(i_subj) = sum(qa_data.ica.final.removed);

        % 9. Offsets
        chan_offsets(i_subj) = mean(abs(qa_data.badchaninfo.offsets.offsets_mV(1:128, :)), 'all');

        % 10. Covariance Matrix
        cov_matrices(:, :, i_subj) = qa_data.channelcov.cov_hf(1:128, 1:128); % cov_hf / cov_clean

        valid_subjects(i_subj) = true;
    end
end
fprintf('Done!\n');

% Filter out missing data
subjects         = subjects(valid_subjects);
N_trials_start   = N_trials_start(valid_subjects);
N_trials_removed = N_trials_removed(valid_subjects);
N_interp_chan    = N_interp_chan(valid_subjects);
N_interp_trl     = N_interp_trl(valid_subjects);
P_emg_left       = P_emg_left(valid_subjects);
V_eog_left       = V_eog_left(valid_subjects);
V_shift          = V_shift(:, valid_subjects);
auto_status      = auto_status(valid_subjects);
ica_redo         = ica_redo(valid_subjects);
cluster_sizes    = cluster_sizes(valid_subjects);
chan_offsets     = chan_offsets(valid_subjects);
ica_removed      = ica_removed(valid_subjects);
cov_matrices     = cov_matrices(:, :, valid_subjects);
NSUB_valid       = length(subjects);
V_shift_med      = median(V_shift, 1)';

% =========================================================================
% Data Quality Index (DQI) Calculation
% =========================================================================
norm_zero = @(x) log1p(abs(x)) ./ (max(log1p(abs(x))) + eps);
Z_chan    = norm_zero(N_interp_chan);
Z_trl     = norm_zero(N_interp_trl);
Z_emg     = norm_zero(P_emg_left);
Z_cluster = norm_zero(cluster_sizes);
Z_ica     = norm_zero(ica_removed);
Z_offset  = norm_zero(chan_offsets);

v_best    = prctile(V_shift_med, 5);
v_worst   = prctile(V_shift_med, 95);
norm_cont = @(x) max(0, min(1, (x - v_best) ./ (v_worst - v_best + eps)));
Z_volt    = norm_cont(V_shift_med);

Z_auto = auto_status;
Z_auto(isnan(Z_auto)) = 0.5;

% Check Interpolated Trials for variance
valid_trl = ~all(isnan(N_interp_trl)) && any(N_interp_trl > 0);

% Calculate composite index
if valid_trl
    DQI = Z_chan + Z_trl + Z_emg + Z_volt + Z_auto;
else
    DQI = Z_chan + Z_emg + Z_volt + Z_auto;
end

% Sort and prepare matrix
[DQI_sorted, sort_idx] = sort(DQI, 'descend');
subjects_sorted = subjects(sort_idx);

% Matrix_Z_Sorted = [Z_chan(sort_idx), Z_emg(sort_idx), Z_volt(sort_idx), ...
%     auto_status(sort_idx), ica_redo(sort_idx), ...
%     Z_cluster(sort_idx), Z_ica(sort_idx), Z_offset(sort_idx)];
% heatmap_labels = {'Interp Chans', 'Leftover EMG', 'Median Voltage', ...
%     'Automagic Status', 'ICA Redo', ...
%     'Interp Elec Cluster Size', 'Num ICs Removed', 'Mean Chan Offset'};

Matrix_Z_Sorted = [Z_chan(sort_idx), Z_emg(sort_idx), Z_volt(sort_idx), ...
    auto_status(sort_idx), ...
    Z_cluster(sort_idx), Z_ica(sort_idx), Z_offset(sort_idx)];
heatmap_labels = {'Interp Chans', 'Leftover EMG', 'Median Voltage', ...
    'Automagic Status', ...
    'Interp Elec Cluster Size', 'Num ICs Removed', 'Mean Chan Offset'};

if valid_trl
    Matrix_Z_Sorted = [Matrix_Z_Sorted(:, 1), Z_trl(sort_idx), Matrix_Z_Sorted(:, 2:end)];
    heatmap_labels  = [heatmap_labels(1), {'Interp Trials'}, heatmap_labels(2:end)];
end

% =========================================================================
% Visualisation 1: The Traffic Light Dashboard
% =========================================================================
fh1 = figure('Name', 'Quality Assurance Dashboard', 'Color', 'w', 'Position', [100, 100, 1100, 800]);
t = tiledlayout(1, 4, 'TileSpacing', 'compact');

nexttile([1 3]);
h = heatmap(heatmap_labels, subjects_sorted, Matrix_Z_Sorted);
custom_cmap  = [linspace(0,1,50)', linspace(0.8,1,50)', linspace(0.2,0.2,50)'];
custom_cmap2 = [linspace(1,0.8,50)', linspace(1,0,50)', linspace(0.2,0.2,50)'];
my_cmap      = [custom_cmap; custom_cmap2];

h.Title       = 'Participant Quality Matrix (Normalised)';
h.Colormap    = my_cmap;
h.ColorLimits = [0 1];
h.XLabel      = 'Quality Metrics (0 = Good, 1 = Bad)';
h.YLabel      = 'Participants (Worst to Best)';

% Tile 4: Composite DQI Bar
nexttile;
barh(DQI_sorted, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
ylim([0.5, NSUB_valid + 0.5]);
yticks([]);
title('Data Quality Index');
xlabel('Total Deviation Score');
warning_threshold = prctile(DQI, 90);
xline(warning_threshold, 'r--', '90th Percentile', 'LabelVerticalAlignment', 'bottom');

% Save Dashboard
plotX = 30; plotY = max(15, NSUB_valid * 0.5);
set(fh1, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh1, fullfile(reports_dir, ['Dashboard_TrafficLight_' myPaths.group]), '-dtiff', '-r300');

% =========================================================================
% Visualisation 2: Distribution Swarmcharts
% =========================================================================
metrics = [{N_trials_start, N_trials_removed}, ...
    {N_interp_chan, P_emg_left, V_eog_left, V_shift_med, cluster_sizes, ica_removed, chan_offsets}];

titles  = [{'Number of Trials (Left)', 'Number of Trials (Removed)'}, ...
    {'Interpolated Channels (%)', 'Leftover EMG (%)', 'Leftover EOG Peak (uV)', 'Median Voltage Shift (uV)', ...
    'Max Bad Elec Cluster Size', 'Number of ICs Removed', 'Mean Channel Offset (mV)'}];

% Dynamically insert Interpolated Trials at position 4 if valid
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

    % KDE Curve
    if std(data) > 1e-6
        [f, x_eval] = ksdensity(data, extraArgs{:});
        f_scaled = (f / max(f)) * 0.5;
        fill(x_eval, f_scaled + 1.3, [0.2 0.4 0.7], 'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.4 0.7]);
    end

    y_base = 1.0;
    y_jitter = y_base + (rand(length(data), 1) - 0.5) * 0.3;
    scatter(data, y_jitter, 25, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.3);

    % Outlier Marking:
    % - For 'Number of Trials (Left)', fewest remaining (left tail) is worst -> mink
    % - For 'Number of Trials (Removed)' and all artifact metrics, highest count (right tail) is worst -> maxk
    num_worst = min(3, length(data));
    if is_trials_left
        [~, worst_idx] = mink(data, num_worst);
    else
        [~, worst_idx] = maxk(data, num_worst);
    end

    scatter(data(worst_idx), y_jitter(worst_idx), 45, 'r', 'filled', 'MarkerEdgeColor', 'k');
    for w = 1:num_worst
        text(data(worst_idx(w)), y_jitter(worst_idx(w)) + 0.18, subs_valid{worst_idx(w)}, ...
            'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
    end

    q = prctile(data, [25 50 75]);
    line([q(2) q(2)], [0.7 1.3], 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5);
    line([q(1) q(3)], [y_base y_base], 'Color', 'k', 'LineWidth', 1.5);
    title(titles{i_metric}, 'FontSize', 12, 'FontWeight', 'bold');

    % Adjust X-limits
    % if contains(titles{i_metric}, '(%)') && ~contains(titles{i_metric}, 'Leftover EMG')
    %     xlim([-0.05 0.5]);
    % else
    %     x_pad = max(std(data) * 0.2, 1e-3);
    %     xlim([min(data) - x_pad, max(data) + x_pad]);
    % end
    x_pad = max(std(data) * 0.2, 1e-3);
    xlim([min(data) - x_pad, max(data) + x_pad]);

    set(gca, 'YTick', [], 'YColor', 'none', 'Box', 'off', 'TickDir', 'out');
    grid on; ax = gca; ax.GridAlpha = 0.1;
end

% Save Distributions
plotX = Ncol * 10; plotY = max(10, Nrow * 8);
set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(reports_dir, ['Dashboard_Distributions_' myPaths.group]), '-dtiff', '-r300');

% =========================================================================
% Visualisation 3: Automated Audit of Voltage Outliers
% =========================================================================
[~, worst_volt_idx] = maxk(V_shift_med, 3);
worst_volt_subjects = subjects(worst_volt_idx);
fprintf('\nTriggering visual audit for the 3 participants with the highest voltage shifts...\n');
audit_voltage_offenders(myPaths, worst_volt_subjects);

% =========================================================================
% Visualisation 4: Channel covariance
% =========================================================================
% CorrelationMatrices2 = CorrelationMatrices;
% load('C:\DATA\MATLAB\myCodes\Preprocessing\files\noisyCov.mat', 'noisyCov');
% CorrelationMatrices2(:, :, 1) = noisyCov; % TEST!
[deviant_indices, fh] = check_channelcov(cov_matrices, subjects);

plotX = 35; plotY = 45;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh, fullfile(reports_dir, ['Summary3_' myPaths.group '_T' num2str(myPaths.visit) '_' myPaths.task  '_' myPaths.proctime]), '-dtiff', '-r400');


fprintf('Dashboards and audits saved successfully to %s\n', reports_dir);

end


function audit_voltage_offenders(myPaths, flagged_subjects)
% Loads the full preprocessed data for flagged subjects and plots a hybrid audit:
% 1. Topoplot of the voltage shift (Where is the noise?)
% 2. PSD of the worst channels (What is the frequency of the noise?)
% 3. Raw trace (What does it look like?)

NSUB = length(flagged_subjects);
reports_dir = fullfile(myPaths.preproc, 'reports');

% 1. Get the dimensions of your primary monitor
% screen_size is [left, bottom, width, height]
screen_size = get(0, 'ScreenSize');
monitor_h = screen_size(4);

% 2. Calculate the desired height, but cap it at 85% of the monitor height
% This ensures the window never "teleports" off the top of the screen.
requested_h = 350 * NSUB;
final_h = min(requested_h, monitor_h * 0.85);

% 3. Calculate a smart 'bottom' position so the window stays visible
% This places the window 100 pixels from the bottom of the screen.
pos_bottom = 100;

% 4. Spawn the figure
fh3 = figure('Name', 'Audit: Spatial & Spectral Check', ...
    'Color', 'w', ...
    'Position', [100, pos_bottom, 1200, final_h]);

for i = 1:NSUB
    subject_id = flagged_subjects{i};
    fprintf('  Auditing %s...\n', subject_id);

    % Define the path to the fully preprocessed continuous data
    path_data = fullfile(myPaths.preproc, subject_id, myPaths.codever, 'data');
    path_data_full = fullfile(path_data, [subject_id '_T' num2str(myPaths.visit) '_' myPaths.task '_cleandata_b.mat']);

    if exist(path_data_full, 'file') == 2
        % Load the heavy EEG data
        temp = load(path_data_full, 'EEG');
        EEG = temp.EEG;

        % Extract the channel-by-channel median voltage shift from your QA data
        volt_shift = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        % Automatically find the 4 worst electrodes driving the high voltage
        [~, worst_idx] = maxk(volt_shift, 4);
        worst_labels = {EEG.chanlocs(worst_idx).labels};
        worst_data = EEG.data(worst_idx, :);

        % -------------------------------------------------------------
        % Plot 1: Spatial Topoplot (Your Method)
        % -------------------------------------------------------------
        nexttile;
        myCmap = brewermap(128, 'Reds');

        % Plot the voltage shift and circle the worst electrodes in black
        topoplot(volt_shift, EEG.chanlocs(1:128), 'maplimits', [prctile(volt_shift, 5), max(volt_shift)], ...
            'headrad', 0.5, 'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'off', ...
            'style', 'map', 'shading', 'interp', ...
            'emarker2', {worst_idx, 'o', 'k', 6, 1});

        title(sprintf('%s: Median Voltage Shift', subject_id));
        hcb = colorbar; hcb.Title.String = '\muV';

        % -------------------------------------------------------------
        % Plot 2: Power Spectral Density (PSD)
        % -------------------------------------------------------------
        nexttile; hold on;

        % Calculate Welch's PSD strictly on the worst channels
        window = EEG.srate * 2;
        noverlap = 0;
        [pxx, f] = pwelch(worst_data', window, noverlap, window, EEG.srate);

        % Plot the mean power across those worst channels
        plot(f, 10*log10(mean(pxx, 2)), 'k', 'LineWidth', 2);

        % Highlight the Alpha Band (8-12 Hz) in light blue
        patch([7 12 12 7], [min(ylim) min(ylim) max(ylim) max(ylim)], [0 0.4 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        title(['PSD of Worst Chans: ' strjoin(worst_labels, ', ')]);
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        xlim([1 70]);
        ylim([-30 30]);
        grid on;

        % -------------------------------------------------------------
        % Plot 3: 5-Second Continuous Trace
        % -------------------------------------------------------------
        nexttile; hold on;

        % Extract a 5-second snippet from the middle of the recording
        mid_point = floor(size(worst_data, 2) / 2);
        win_samples = 5 * EEG.srate;
        start_idx = max(1, mid_point - floor(win_samples/2));
        end_idx = min(size(worst_data, 2), start_idx + win_samples - 1);

        t_vec = (0:(end_idx - start_idx)) / EEG.srate;
        plot(t_vec, worst_data(:, start_idx:end_idx)', 'LineWidth', 1.2);

        title('5s Raw Trace of Worst Chans');
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        ylim([-100 100]);
        grid on;

    else
        % Handle missing files smoothly
        nexttile([1 3]);
        text(0.5, 0.5, sprintf('Data missing for %s', subject_id), 'HorizontalAlignment', 'center', 'FontSize', 12);
        axis off;
    end
end

% Save the Audit Figure
plotX = 40; plotY = max(10, NSUB * 8);
set(fh3, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh3, fullfile(reports_dir, ['Dashboard_Audit_Voltage_' myPaths.group]), '-dtiff', '-r300');

end