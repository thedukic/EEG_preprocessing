function report_final_new(myPaths, reports_dir)
% =========================================================================
% Script for reporting on EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
% =========================================================================

fprintf('\n==================================================================\n');
fprintf('%s: Generating the final quality dashboard.\n', myPaths.group);
fprintf('==================================================================\n');

% -------------------------------------------------------------------------
% Path Standardisation
% -------------------------------------------------------------------------
% Ensure preproc paths are stored as a cell array for consistent handling
if ischar(myPaths.preproc) || isstring(myPaths.preproc)
    preproc_paths = {char(myPaths.preproc)};
elseif iscell(myPaths.preproc)
    preproc_paths = myPaths.preproc;
else
    error('myPaths.preproc must be a character array, string, or cell array.');
end

% Set reports directory to the first path in the list
% reports_dir = fullfile(preproc_paths{1}, 'reports');
if exist(reports_dir,'dir') ~= 7, mkdir(reports_dir); end

% -------------------------------------------------------------------------
% Pre-scan Folders for Allocation
% -------------------------------------------------------------------------
total_sub_count = 0;
path_subject_map = cell(length(preproc_paths), 2);

for p_idx = 1:length(preproc_paths)
    this_path = preproc_paths{p_idx};
    these_subjects = list_participants(this_path, {});
    
    if ~isempty(these_subjects)
        if isrow(these_subjects), these_subjects = these_subjects'; end
        path_subject_map{p_idx, 1} = this_path;
        path_subject_map{p_idx, 2} = these_subjects;
        total_sub_count = total_sub_count + length(these_subjects);
    end
end

if total_sub_count == 0
    fprintf('No subjects found in the provided directories. Exiting.\n');
    return;
end

% -------------------------------------------------------------------------
% Strict Folder-by-Folder Processing
% -------------------------------------------------------------------------
% Preallocate arrays based on total possible subjects
N_interp_chan  = NaN(total_sub_count, 1);
N_interp_trl   = NaN(total_sub_count, 1);
N_emg_left     = NaN(total_sub_count, 1);
V_shift        = NaN(128, total_sub_count);
Gamma_CV       = NaN(total_sub_count, 1); 
Auto_Status    = NaN(total_sub_count, 1); 

% Tracking arrays for the final dashboard
subjects_master = cell(total_sub_count, 1);
paths_master    = cell(total_sub_count, 1);

fprintf('Loading QA metrics from target folders... \n');

global_idx = 0; % Counter for successfully loaded files

for p_idx = 1:size(path_subject_map, 1)
    current_path = path_subject_map{p_idx, 1};
    current_subjects = path_subject_map{p_idx, 2};
    
    if isempty(current_path), continue; end
    
    % Set tempPaths to ONLY look in the current folder
    tempPaths = myPaths;
    tempPaths.preproc = current_path;
    
    for i_subj = 1:length(current_subjects)
        subject_id = current_subjects{i_subj};
        
        % preproc_folders_subject is strictly fed the current path
        subject = preproc_folders_subject(subject_id, tempPaths, 2);
        qa_name = fullfile(subject.qa, subject.qametrics{1});

        if exist(qa_name, 'file')
            global_idx = global_idx + 1;
            
            % Record the subject and its explicit source path
            subjects_master{global_idx} = subject_id;
            paths_master{global_idx} = current_path;
            
            % Load Data
            load(qa_name, 'qa_data');
            qa_data = report_issues(qa_data);

            % 1. Interpolated Channels (%)
            bad_chans = qa_data.badchaninfo.badElectrodes;
            N_interp_chan(global_idx) = length(bad_chans) / 128;

            % 2. Interpolated Trials (%)
            if isfield(qa_data.epochRejections, 'InterpTrialInfo') && ~isempty(qa_data.epochRejections.InterpTrialInfo)
                total_trials = qa_data.issues_to_check.NumberTrials1;
                interp_epochs = qa_data.epochRejections.interpEpochs;
                N_interp_trl(global_idx) = interp_epochs / total_trials;
            else
                N_interp_trl(global_idx) = 0;
            end

            % 3. EMG Leftover
            if isfield(qa_data.leftovers, 'muscle2')
                N_emg_left(global_idx) = qa_data.leftovers.muscle2;
            end

            % 4. Voltage Shift
            V_shift(:, global_idx) = qa_data.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

            % 5. Gamma Spread
            if isfield(qa_data, 'psd_gamma_cv')
                Gamma_CV(global_idx) = qa_data.psd_gamma_cv;
            end

            % 6. Extract Automagic Status
            status_str = lower(strtrim(qa_data.automagicmetrics.status));
            switch status_str
                case 'good', Auto_Status(global_idx) = 0;   % Green
                case 'ok',   Auto_Status(global_idx) = 0.5; % Yellow
                case 'bad',  Auto_Status(global_idx) = 1;   % Red
                otherwise,   Auto_Status(global_idx) = NaN;
            end
        end
    end
end
fprintf('Done! Successfully loaded %d datasets.\n', global_idx);

if global_idx == 0
    fprintf('No valid QA files found. Exiting.\n');
    return;
end

% -------------------------------------------------------------------------
% Trim arrays down to the successful loads
% -------------------------------------------------------------------------
subjects      = subjects_master(1:global_idx);
subject_paths = paths_master(1:global_idx);
N_interp_chan = N_interp_chan(1:global_idx);
N_interp_trl  = N_interp_trl(1:global_idx);
N_emg_left    = N_emg_left(1:global_idx);
V_shift       = V_shift(:, 1:global_idx);
Auto_Status   = Auto_Status(1:global_idx);
NSUB_valid    = global_idx;

V_shift_med = median(V_shift, 1)';

% =========================================================================
% Data Quality Index (DQI) Calculation
% =========================================================================
norm_zero = @(x) log1p(abs(x)) ./ (max(log1p(abs(x))) + eps);
Z_chan = norm_zero(N_interp_chan);
Z_trl  = norm_zero(N_interp_trl);
Z_emg  = norm_zero(N_emg_left);

v_best = prctile(V_shift_med, 5);
v_worst = prctile(V_shift_med, 95);
norm_cont = @(x) max(0, min(1, (x - v_best) ./ (v_worst - v_best + eps)));
Z_volt = norm_cont(V_shift_med);

Z_auto = Auto_Status;
Z_auto(isnan(Z_auto)) = 0.5;

DQI = Z_chan + Z_trl + Z_emg + Z_volt + Z_auto;

[DQI_sorted, sort_idx] = sort(DQI, 'descend');
subjects_sorted = subjects(sort_idx);
Matrix_Z_Sorted = [Z_chan(sort_idx), Z_trl(sort_idx), Z_emg(sort_idx), Z_volt(sort_idx), Auto_Status(sort_idx)];

% =========================================================================
% Visualisation 1: The Traffic Light Dashboard
% =========================================================================
fh1 = figure('Name', 'Quality Assurance Dashboard', 'Color', 'w', 'Position', [100, 100, 1100, 800]);
t = tiledlayout(1, 4, 'TileSpacing', 'compact');

nexttile([1 3]);
h = heatmap({'Interp Chans', 'Interp Trials', 'Leftover EMG', 'Median Voltage', 'Automagic Status'}, ...
    subjects_sorted, Matrix_Z_Sorted);

custom_cmap = [linspace(0,1,50)', linspace(0.8,1,50)', linspace(0.2,0.2,50)'];
custom_cmap2 = [linspace(1,0.8,50)', linspace(1,0,50)', linspace(0.2,0,50)'];
my_cmap = [custom_cmap; custom_cmap2];

h.Title = 'Participant Quality Matrix (Normalised)';
h.Colormap = my_cmap; 
h.ColorLimits = [0 1]; 
h.XLabel = 'Quality Metrics (0 = Good, 1 = Bad)';
h.YLabel = 'Participants (Worst to Best)';

nexttile;
barh(DQI_sorted, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse'); 
ylim([0.5 NSUB_valid+0.5]);
yticks([]); 
title('Data Quality Index');
xlabel('Total Deviation Score');

warning_threshold = prctile(DQI, 90);
xline(warning_threshold, 'r--', '90th Percentile', 'LabelVerticalAlignment', 'bottom');

plotX = 30; plotY = max(15, NSUB_valid * 0.5); 
set(fh1, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh1, fullfile(reports_dir, ['Dashboard_TrafficLight_' myPaths.group]), '-dtiff', '-r300');

% =========================================================================
% Visualisation 2: Distribution Swarmcharts
% =========================================================================
fh2 = figure('Name', 'Cohort Distributions', 'Color', 'w', 'Position', [150, 150, 1000, 700]);
t2 = tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

metrics = {N_interp_chan, N_interp_trl, N_emg_left, V_shift_med};
titles = {'Interpolated Channels (%)', 'Interpolated Trials (%)', 'Leftover EMG (%)', 'Median Voltage Shift (uV)'};

for i_subj = 1:4
    nexttile; hold on;

    raw_data = metrics{i_subj};
    valid_mask = ~isnan(raw_data);
    data = raw_data(valid_mask);
    subs_valid = subjects(valid_mask); 

    if isempty(data), continue; end

    extraArgs = {};
    if i_subj < 4
        low_s  = -eps;
        high_s = max(data) * 1.2 + eps;
        extraArgs = {'BoundaryCorrection', 'reflection', 'Support', [low_s, high_s]};
    else
        low_s  = min(data) - (std(data)*0.5);
        high_s = max(data) + (std(data)*0.5);
        extraArgs = {'Support', [low_s, high_s]};
    end

    [f, x_eval] = ksdensity(data, extraArgs{:});

    f_scaled = (f / max(f)) * 0.5;
    fill(x_eval, f_scaled + 1.3, [0.2 0.4 0.7], 'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.4 0.7]);

    y_base = 1.0;
    y_jitter = y_base + (rand(length(data), 1) - 0.5) * 0.3;
    scatter(data, y_jitter, 25, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.3);

    [~, worst_idx] = maxk(data, 3); 
    scatter(data(worst_idx), y_jitter(worst_idx), 45, 'r', 'filled', 'MarkerEdgeColor', 'k');

    for w = 1:length(worst_idx)
        text(data(worst_idx(w)), y_jitter(worst_idx(w)) + 0.18, subs_valid{worst_idx(w)}, ...
            'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
    end

    q = prctile(data, [25 50 75]);
    line([q(2) q(2)], [0.7 1.3], 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5); 
    line([q(1) q(3)], [y_base y_base], 'Color', 'k', 'LineWidth', 1.5);    

    title(titles{i_subj}, 'FontSize', 12, 'FontWeight', 'bold');
    if i_subj < 4, xlim([0 0.5]); else, xlim([min(data)-5, max(data)+5]); end

    set(gca, 'YTick', [], 'YColor', 'none', 'Box', 'off', 'TickDir', 'out');
    grid on; ax = gca; ax.GridAlpha = 0.1;
end

plotX = 25; plotY = 15;
set(fh2, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh2, fullfile(reports_dir, ['Dashboard_Distributions_' myPaths.group]), '-dtiff', '-r300');

% =========================================================================
% Visualisation 3: Automated Audit of Voltage Outliers
% =========================================================================
[~, worst_volt_idx] = maxk(V_shift_med, 3);
worst_volt_subjects = subjects(worst_volt_idx);
worst_volt_paths = subject_paths(worst_volt_idx); 

fprintf('\nTriggering visual audit for the 3 participants with the highest voltage shifts...\n');
audit_voltage_offenders(myPaths, worst_volt_subjects, worst_volt_paths, reports_dir);

fprintf('Dashboards and audits saved successfully to %s\n', reports_dir);
end


function audit_voltage_offenders(myPaths, flagged_subjects, flagged_paths, reports_dir)
% Loads the full preprocessed data for flagged subjects using their exact stored path

NSUB = length(flagged_subjects);

screen_size = get(0, 'ScreenSize');
monitor_h = screen_size(4);
requested_h = 350 * NSUB;
final_h = min(requested_h, monitor_h * 0.85);
pos_bottom = 100;

fh3 = figure('Name', 'Audit: Spatial & Spectral Check', ...
    'Color', 'w', ...
    'Position', [100, pos_bottom, 1200, final_h]);

for i = 1:NSUB
    subject_id = flagged_subjects{i};
    subject_path = flagged_paths{i}; 
    fprintf('  Auditing %s from %s...\n', subject_id, subject_path);

    % Go directly to the correct path that was mapped during the initial scan
    path_data = fullfile(subject_path, subject_id, myPaths.codever, 'data');
    path_data_full = fullfile(path_data, [subject_id '_T' num2str(myPaths.visit) '_' myPaths.task '_cleandata_b.mat']);
        
    if exist(path_data_full, 'file') == 2
        temp = load(path_data_full, 'EEG');
        EEG = temp.EEG;

        volt_shift = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        [~, worst_idx] = maxk(volt_shift, 4);
        worst_labels = {EEG.chanlocs(worst_idx).labels};
        worst_data = EEG.data(worst_idx, :);

        % Plot 1
        nexttile;
        myCmap = brewermap(128, 'Reds');

        topoplot(volt_shift, EEG.chanlocs(1:128), 'maplimits', [prctile(volt_shift, 5), max(volt_shift)], ...
            'headrad', 0.5, 'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'off', ...
            'style', 'map', 'shading', 'interp', ...
            'emarker2', {worst_idx, 'o', 'k', 6, 1});

        title(sprintf('%s: Median Voltage Shift', subject_id), 'Interpreter', 'none');
        hcb = colorbar; hcb.Title.String = '\muV';

        % Plot 2
        nexttile; hold on;
        window = EEG.srate * 2;
        noverlap = 0;
        [pxx, f] = pwelch(worst_data', window, noverlap, window, EEG.srate);

        plot(f, 10*log10(mean(pxx, 2)), 'k', 'LineWidth', 2);
        patch([7 12 12 7], [min(ylim) min(ylim) max(ylim) max(ylim)], [0 0.4 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        title(['PSD of Worst Chans: ' strjoin(worst_labels, ', ')], 'Interpreter', 'none');
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        xlim([1 70]);
        ylim([-30 30]);
        grid on;

        % Plot 3
        nexttile; hold on;
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
        nexttile([1 3]);
        text(0.5, 0.5, sprintf('Data missing for %s', subject_id), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'none');
        axis off;
    end
end

plotX = 40; plotY = max(10, NSUB * 8);
set(fh3, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh3, fullfile(reports_dir, ['Dashboard_Audit_Voltage_' myPaths.group]), '-dtiff', '-r300');

end