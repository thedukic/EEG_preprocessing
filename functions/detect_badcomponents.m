function EEG = detect_badcomponents(EEG, EXT, EMG, cfg)

threshold_ext = 3;
% threshold_ext_blink    = 3;
% threshold_ext_saccade  = 3;
% threshold_ext_heart    = 3;

threshold_temp_blink   = 0.85;
threshold_temp_saccade = 0.85;
threshold_temp_heart   = 0.85;

threshold_ctps_pk = 20; % paper: 20

threshold_power_int = 10; % empirical: 10

% threshold_channel_xxx = [3 2000]; % empirical: [3 2000]
threshold_channel_kurt = 20; % empirical: 20

%% ========================================================================
fprintf('\n================================\n');
fprintf('Detecting bad ICs\n');
fprintf('================================\n');

% Double-check
EEG = eeg_checkset(EEG, 'ica');
EEG.icaact = [];

% Make sure IC activations are present
data_ica   = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :);
num_ica    = size(data_ica,1);

% Smoothen the ICs to remove 'noise' and improve the correlation
icawinv = EEG.icawinv;
% icawinvSmooth = estimate_invlaplacian(icawinv, EEG.chanlocs, 1);

% Load individual/template ICA templates
templates_ica = load_ictemplateweights(EEG);

%% ========================================================================
% 1. ICLabel
fprintf('\n--------------------------------\n');
fprintf('ICLabel\n');
fprintf('--------------------------------\n');

% Build active matrix on the fly
active_thresholds = cfg.ica.iclabel.base;

switch lower(EEG.ALSUTRECHT.subject.task)
    case {'mmn', 'sart'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.erp;
    case {'mt'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.mt;
    case {'rs', 'eo', 'ec'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.rs;
end

% Detect bad ICs
EEG = iclabel(EEG);
EEG = pop_icflag(EEG, active_thresholds);

% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel.bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel.clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

%% ========================================================================
% 2. Correlation with external channels (ECG / VEOG / HEOG)
fprintf('\n--------------------------------\n');
fprintf('ECG/VEOG/HEOG: Correlations with the external electrodes\n');
fprintf('--------------------------------\n');

channel_ecg  = find(strcmp({EXT.chanlocs.labels}, 'ECG'));
channel_veog = find(strcmp({EXT.chanlocs.labels}, 'VEOG'));
channel_heog = find(strcmp({EXT.chanlocs.labels}, 'HEOG'));

% Check if ECG signal is recorded
if isempty(channel_ecg)
    fprintf('Warning: ECG signal was not recorded, so some heart IC detections cannot be done.\n');
    data_ext_ecg = [];
else
    data_ext_ecg = EXT.data(channel_ecg, :);
end

% EOG signlas
data_ext_eog = EXT.data([channel_veog channel_heog], :);

% Temporarily filter for better detection
% https://mne.tools/stable/generated/mne.preprocessing.create_eog_epochs.html
% https://ieeexplore.ieee.org/document/4536072
% VEOG: 0-10 Hz
% HEOG: 0-30 Hz
% ECG:  8-16 Hz / 10-20 Hz

% EOG
% Keep very low freq
[bh_eog, ah_eog] = butter(2, 0.3/(EEG.srate/2), 'high');
[bl_eog, al_eog] = butter(2, 10/(EEG.srate/2), 'low');
% [bh_eog, ah_eog] = butter(2, 1/(EEG.srate/2), 'high');
% [bl_eog, al_eog] = butter(2, 20/(EEG.srate/2), 'low');

% Filter ICA for EOG correlation
data_ica_eog = do_filteringcore(bl_eog, al_eog, data_ica, EEG.event, EEG.srate);
data_ica_eog = do_filteringcore(bh_eog, ah_eog, data_ica_eog, EEG.event, EEG.srate);
data_ica_eog = data_ica_eog';

% Filter EXT for EOG correlation
data_ext_eog = do_filteringcore(bl_eog, al_eog, data_ext_eog, EEG.event, EEG.srate);
data_ext_eog = do_filteringcore(bh_eog, ah_eog, data_ext_eog, EEG.event, EEG.srate);

% ECG
if ~isempty(data_ext_ecg)
    [bh_ecg, ah_ecg] = butter(2, 10/(EEG.srate/2), 'high');
    [bl_ecg, al_ecg] = butter(2, 20/(EEG.srate/2), 'low');

    % Filter ICA for ECG correlation
    data_ica_ecg = do_filteringcore(bl_ecg, al_ecg, data_ica, EEG.event, EEG.srate);
    data_ica_ecg = do_filteringcore(bh_ecg, ah_ecg, data_ica_ecg, EEG.event, EEG.srate);
    data_ica_ecg = data_ica_ecg';

    % Filter EXT for ECG correlation
    data_ext_ecg = do_filteringcore(bl_ecg, al_ecg, data_ext_ecg, EEG.event, EEG.srate);
    data_ext_ecg = do_filteringcore(bh_ecg, ah_ecg, data_ext_ecg, EEG.event, EEG.srate);
else
    data_ica_ecg = NaN * data_ica';
    data_ext_ecg = NaN * data_ext_eog(1, :);
end

% Combine: ECG / VEOG / HEOG
% threshold_ext = [threshold_ext_heart, threshold_ext_blink, threshold_ext_saccade];
labels_ext    = {'ECG', 'VEOG', 'HEOG'};
data_ext_all  = [data_ext_ecg; data_ext_eog]';
clearvars data_ext_ecg data_ext_eog

% Correlation: ECG
corr_ecg = corr(data_ica_ecg.^2, abs(data_ext_all(:, 1)), "type", "Spearman");
% corr_ecg = corr(data_ica_ecg, data_ext_all(:, 1), "type", "Spearman");

% Correlation: VEOG
corr_veog = corr(data_ica_eog, data_ext_all(:, 2), "type", "Spearman");

% Correlation: HEOG
corr_heog = corr(data_ica_eog, data_ext_all(:, 3), "type", "Spearman");

% Combine: ECG / VEOG / HEOG
% corr_ext = abs([corr_ecg, corr_veog, corr_heog]);
corr_ext = abs(zscore([corr_ecg, corr_veog, corr_heog]));

num_ext = size(corr_ext, 2);
assert(length(labels_ext) == num_ext);
done_ext = false(1, num_ext);

for i_ext = 1:num_ext
    switch i_ext
        case 1
            % ECG, only 1 IC
            corr_tmp = corr_ext(:, i_ext);
            if ~any(isnan(corr_tmp))
                done_ext(i_ext) = true;

                bad_ic_ecg = find(corr_tmp > threshold_ext);

                if length(bad_ic_ecg) > 1
                    [~, bad_ic_mostlikely] = max(corr_tmp(bad_ic_ecg));
                    bad_ic_ecg = bad_ic_ecg(bad_ic_mostlikely);
                end
            else
                bad_ic_ecg = [];
            end

        case 2
            % VEOG, allow multiple
            done_ext(i_ext) = true;

            corr_tmp = corr_ext(:, i_ext);
            bad_ic_veog = find(corr_tmp > threshold_ext);

        case 3
            % HEOG, only 1 IC with max corr
            done_ext(i_ext) = true;

            corr_tmp = corr_ext(:, i_ext);
            bad_ic_heog = find(corr_tmp > threshold_ext);

            if length(bad_ic_heog) > 1
                [maxValue, bad_ic_mostlikely] = max(corr_tmp(bad_ic_heog));
                bad_ic_heog = bad_ic_heog(bad_ic_mostlikely);
            end
    end

    fprintf('%s: max{abs(zscore(R))} = %1.2f\n', labels_ext{i_ext}, max(corr_tmp));
end

% Combine: ECG / VEOG / HEOG
bad_ic = [bad_ic_ecg(:); bad_ic_veog(:); bad_ic_heog(:)];
bad_ic_type = [ones(length(bad_ic_ecg), 1); 2*ones(length(bad_ic_veog), 1); 3*ones(length(bad_ic_heog), 1)];

% Report
fprintf('\n');
for i_ext = 1:num_ext
    if done_ext(i_ext)
        if any(bad_ic_type == i_ext)
            fprintf('%s ICs (N = %d) were identified using EXT channel correlation (Z > %1.1f).\n', labels_ext{i_ext}, sum(bad_ic_type == i_ext), threshold_ext);
            str = strjoin(string(bad_ic(bad_ic_type == i_ext)), ', ');
            fprintf('%s ICs: %s\n',labels_ext{i_ext}, str);
        else
            fprintf('No %s ICs were identified using EXT channel correlation (Z > %1.1f).\n', labels_ext{i_ext}, threshold_ext);
        end
    else
        fprintf('No %s ICs were checked using EXT channel correlation.\n', labels_ext{i_ext});
    end
end

% Log
EEG.ALSUTRECHT.ica.corr.corr = single(corr_ext);
EEG.ALSUTRECHT.ica.corr.bics = bad_ic;
EEG.ALSUTRECHT.ica.corr.cvec = bad_ic_type;
EEG.ALSUTRECHT.ica.corr.clss = labels_ext;

% -------------------------------------
% Plot
% -------------------------------------
% % Plot (Barplot approach)
% fh = figure('Color', 'w', 'Position', [100, 100, 850, 450], 'Visible', cfg.figure.visible);
%
% % Create a 2x1 tiled chart layout
% % 'TileSpacing', 'tight' dramatically reduces the vertical gap between plots
% % 'Padding', 'compact' reduces the outer white margins around the figure
% t = tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
% ax1 = nexttile(t);
%
% % Create a grouped bar chart for discrete component values
% z_scores = abs(zscore(corr_ext));
% b = bar(1:num_ics, z_scores, 'grouped', 'EdgeColor', 'none');
%
% % Apply a refined, earthy color palette
% b(1).FaceColor = [0.25, 0.40, 0.50]; % Muted Slate Blue (ECG)
% b(2).FaceColor = [0.30, 0.50, 0.35]; % Olive Green (VEOG)
% b(3).FaceColor = [0.55, 0.35, 0.55]; % Muted Purple (HEOG)
%
% hold on;
%
% % Add statistical threshold markers (e.g., z-score = 2.5)
% yline(threshold_ext, '--', 'Color', [0.7, 0.2, 0.2], 'LineWidth', 1.2);
%
% % Style the plot axes
% grid on;
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Layer', 'top');
% set(gca, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);
%
% % Labels and boundaries
% xlabel('Independent Components (ICs)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Correlation Strength (Z-Score)', 'FontSize', 12, 'FontWeight', 'bold');
% title('ICA Component Correlation with External Channels', 'FontSize', 14, 'FontWeight', 'bold');
% xlim([0.5, num_ics + 0.5]);
% xticks(1:num_ics);
%
% % Clean legend configuration
% legend({'ECG', 'VEOG', 'HEOG', 'Threshold'}, 'Location', 'northeast');
% hold off;

% Plot (Heatmap approach)
fh = figure('Color', 'w', 'Position', [100, 100, 850, 450], 'Visible', cfg.figure.visible);

t = tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
ax1 = nexttile(t);

% Calculate absolute Z-scores matrix (num_ext x num_ics)
z_scores = abs(zscore(corr_ext));

% Display heatmap (transpose so num_ics is on X-axis and num_ext on Y-axis)
imagesc(1:num_ica, 1:size(z_scores, 2), z_scores');

% Apply Blue-Red colormap (Blue = low z-score, Red = high z-score)
colormap(ax1, brewermap([], 'Reds'));
c = colorbar;
c.Label.String = 'Correlation Strength (|Z-Score|)';

% Set axis properties
set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11, 'YDir', 'normal');
xlabel('Independent Components (ICs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('External Channels', 'FontSize', 12, 'FontWeight', 'bold');
title('ICA Component Correlation with External Channels', 'FontSize', 14, 'FontWeight', 'bold');

% Ticks setup
xticks(1:num_ica);
yticks(1:3);
yticklabels({'ECG', 'VEOG', 'HEOG'});

plotX = max(20, num_ica * 0.8 + 5);
plotY = plotX / 1.6;
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_external'], [plotX plotY]);

%% ========================================================================
% 3. Cross-trial phase statistics for ECG detection
fprintf('\n--------------------------------\n');
fprintf('ECG: Cross-trial phase statistics\n');
fprintf('--------------------------------\n');

% Initialise
pulse_estimate = NaN;
V      = NaN(num_ica, 1);
pK     = NaN(num_ica, 1);
bad_ic = [];
is_ecg = NaN;
stats  = NaN;

if ~isempty(channel_ecg)
    % Detect ECG
    ecg_window_ms = [-200 250];
    [ecg_mask, ecg_epoch, ~, ~, ~, pulse_estimate] = detect_ecg(EXT, ecg_window_ms, cfg.figure.visible);

    if ~isnan(ecg_mask)
        % #1 CTPS
        [V, pK] = my_ctps(data_ica, ecg_epoch, EEG.event, EEG.srate);
        bad_ic = find(pK >= threshold_ctps_pk);

        if ~isempty(bad_ic)
            fprintf('ECG ICs (N = %d) were identified using cross-trial phase statistics (pK >= %1.1f).\n', length(bad_ic), threshold_ctps_pk);
            str = strjoin(string(bad_ic), ', ');
            fprintf('ECG ICs: %s\n', str);
        else
            fprintf('No ECG ICs were identified using cross-trial phase statistics (pK >= %1.1f).\n', threshold_ctps_pk);
        end

        % #2 ERP
        % Configure with epoch window and threshold
        cfg_ecg = [];
        cfg_ecg.epoch_window_ms = ecg_window_ms;
        cfg_ecg.min_corr        = 0.8;
        cfg_ecg.min_snr         = 3.5;
        cfg_ecg.max_lag_ms      = 30;
        cfg_ecg.ecg_chan        = 'ECG';
        cfg_ecg.do_plot         = true;
        [is_ecg, stats, fh] = check_ic_ecg_erp(EEG, EXT, ecg_epoch, cfg_ecg);
        save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_ecg_erp'], [20 15]);

    else
        fprintf('Skipping. The ECG signal is corrupted or too noisy.\n');
    end
else
    fprintf('Skipping. The ECG signal is not recorded.\n');
end

% Log
EEG.ALSUTRECHT.ica.heart.ctps.v     = V;
EEG.ALSUTRECHT.ica.heart.ctps.pk    = pK;
EEG.ALSUTRECHT.ica.heart.ctps.bics  = bad_ic;
EEG.ALSUTRECHT.ica.heart.ctps.cvec  = ones(length(EEG.ALSUTRECHT.ica.heart.ctps.bics), 1);
EEG.ALSUTRECHT.ica.heart.ctps.clss  = 'ECG';
EEG.ALSUTRECHT.ica.heart.erp.is_ecg = is_ecg;
EEG.ALSUTRECHT.ica.heart.erp.stats  = stats;

EEG.ALSUTRECHT.subject.heartrate = pulse_estimate;

% -------------------------------------
% Plot
% -------------------------------------
fh = figure('Color', 'w', 'Position', [100, 100, 850, 600], 'Visible', cfg.figure.visible);

% Create a 2x1 tiled chart layout
% 'TileSpacing', 'tight' dramatically reduces the vertical gap between plots
% 'Padding', 'compact' reduces the outer white margins around the figure
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- TOP PANEL: pK Value ---
ax1 = nexttile(t);
bar(1:num_ica, pK, 'FaceColor', [0.25, 0.40, 0.50], 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on;

% Add the critical paper threshold at pK = 20
line_pk = yline(20, '--', 'Color', [0.7, 0.2, 0.2], 'LineWidth', 1.5);

% Styling top panel
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Layer', 'top');
set(gca, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);
ylabel('Significance Score (p_K)', 'FontSize', 12, 'FontWeight', 'bold');
title('ICA Component Phase Locking to ECG R-Peaks (CTPS)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0.5, num_ica + 0.5]);
xticks(1:num_ica);
% xticklabels({}); % Hide x-labels on top to save space

% Add a clean legend for the threshold
legend(line_pk, 'Artifact Threshold (p_K \geq 20)', 'Location', 'NorthEast', 'Box', 'off');

% --- BOTTOM PANEL: Kuiper Index V ---
ax2 = nexttile(t);
bar(1:num_ica, V, 'FaceColor', [0.45, 0.45, 0.45], 'EdgeColor', 'none', 'BarWidth', 0.6);
grid on;

% Styling bottom panel
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Layer', 'top');
set(gca, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);
xlabel('Independent Components (ICs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Kuiper Index (V)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5, num_ica + 0.5]);
ylim([0, 1]); % V is naturally bounded between 0 and 1
xticks(1:num_ica);

% Link the x-axes so zooming or panning shifts both panels together
linkaxes([ax1, ax2], 'x');

plotX = max(20, num_ica * 0.8 + 5);
plotY = plotX / 1.6;
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_ctps'], [plotX plotY]);

%% ========================================================================
% 4. Detect EMG ICs using freq slopes
fprintf('\n--------------------------------\n');
fprintf('EMG: Power slopes\n');
fprintf('--------------------------------\n');

options = [];
options.Freq_to_compute       = [1 100];
options.muscleFreqEx          = 50 + 2*[-1 1]; % Line freq +-bandwith
options.muscleFreq1           = cfg.emg.slope_freq_1; % e.g., [7 45]
options.muscleFreq2           = cfg.emg.slope_freq_2; % e.g., [40 70]
options.muscleSlopeThreshold1 = cfg.emg.slope_threshold_1;
options.muscleSlopeThreshold2 = cfg.emg.slope_threshold_2;

% Calculate pwelch to enable detection of log-freq log-power slopes, indicative of muscle activity
[pow, frefAll] = pwelch(data_ica', size(data_ica, 2), [], size(data_ica, 2), EEG.srate);

pow = pow';
frefAll = frefAll';

% Calculate FFT bins
freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
fftBins = zeros(size(pow, 1), size(freq, 2));
for i_freq = 1:length(freq)
    [~, index1] = min(abs(frefAll-((freq(1,i_freq)-0.25))));
    [~, index2] = min(abs(frefAll-((freq(1,i_freq)+0.25))));
    fftBins(:, i_freq) = mean(pow(:, index1:index2), 2); % creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
end

% figure;
% nexttile; plot(fp, pxx(:,1));
% nexttile; plot(freq, fftBins(1,:));

% Better muscle comp_number identification using Dual-Window Slopes:
slope_muscle = NaN(num_ica, 2);

for i_ic = 1:num_ica
    % ---------------------------------------------------------------------
    % 1. Extract data for the Broad Window (e.g., 7-45 Hz)
    % ---------------------------------------------------------------------
    if ~isempty(options.muscleFreq1)
        [~, fin1] = min(abs(options.muscleFreq1(1) - freq));
        [~, fin2] = min(abs(options.muscleFreq1(2) - freq));
        freqHz_broad = freq(fin1:fin2);
        freqPow_broad = fftBins(i_ic, fin1:fin2);
    else
        freqHz_broad = freq;
        freqPow_broad = fftBins(i_ic, :);
    end

    % Exclude line noise from broad window
    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_broad));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_broad));
        if fex1 <= length(freqHz_broad) % Only remove if line noise is within window
            freqHz_broad(fex1:fex2) = [];
            freqPow_broad(fex1:fex2) = [];
        end
    end

    % Fit linear regression to broad window
    p_broad = polyfit(log10(freqHz_broad), log10(freqPow_broad), 1);
    slope_broad = p_broad(1);

    % ---------------------------------------------------------------------
    % 2. Extract data for the High-Frequency Window (e.g., 30-90 Hz)
    % ---------------------------------------------------------------------
    [~, fhi1] = min(abs(options.muscleFreq2(1) - freq));
    [~, fhi2] = min(abs(options.muscleFreq2(2) - freq));
    freqHz_high = freq(fhi1:fhi2);
    freqPow_high = fftBins(i_ic, fhi1:fhi2);

    % Exclude line noise from high window
    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_high));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_high));
        if fex1 <= length(freqHz_high) % Ensure line noise is actually inside this range
            freqHz_high(fex1:fex2) = [];
            freqPow_high(fex1:fex2) = [];
        end
    end

    % Fit linear regression to high-frequency window
    p_high = polyfit(log10(freqHz_high), log10(freqPow_high), 1);
    slope_high = p_high(1);

    % ---------------------------------------------------------------------
    % 3. Combine Logic: Take the most severe (highest) slope
    % ---------------------------------------------------------------------
    % If either the broad slope OR the high-frequency slope is heavily positive,
    % the maximum value will catch it.
    slope_muscle(i_ic, :) = [slope_broad, slope_high];
end

% Detect EMG components
badIC1 = find(slope_muscle(:, 1) > options.muscleSlopeThreshold1 | slope_muscle(:, 2) > options.muscleSlopeThreshold2);

% Combine
bad_ic = badIC1(:);

% Report
report_detections('Muscle', 'power slope', 'slope', bad_ic, options.muscleSlopeThreshold1);

% Log
EEG.ALSUTRECHT.ica.muscle.slope.slope = slope_muscle(:);
EEG.ALSUTRECHT.ica.muscle.slope.bics  = bad_ic(:);
EEG.ALSUTRECHT.ica.muscle.slope.cvec  = ones(sum(EEG.ALSUTRECHT.ica.muscle.slope.bics),1);
EEG.ALSUTRECHT.ica.muscle.slope.clss  = 'EMG';

% % -------------------------------------
% % Plot
% % -------------------------------------
% fh = figure('Color', 'w', 'Position', [100, 100, 950, 650], 'Visible', cfg.figure.visible);
%
% % Create a 2x2 grid with tight spacing and compact padding to eliminate white space
% t = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
%
% % Define consistent style parameters
% bar_color = [0.45, 0.45, 0.45]; % Professional neutral gray for baseline metrics
% ratio_color = [0.55, 0.20, 0.20]; % Distinct dark red for the final selection ratio
% font_name = 'Helvetica';
%
% % --- TILE 1: Negative Integral (Z-score) ---
% nexttile(t);
% data_plot = zscore(neg_integral);
% bar(1:num_ics, data_plot, 'FaceColor', bar_color, 'EdgeColor', 'none', 'BarWidth', 0.6);
% grid on;
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Box', 'off', 'FontName', font_name, 'FontSize', 10);
% title('Spectral Negative Integral (Z-Score)', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('Z-Score');
% xlim([0.5, num_ics + 0.5]);
% xticks(1:num_ics);
% xticklabels({}); % Hide x-labels to keep the top row clean
% ylim_max = max(4, round(max(abs(data_plot))));
% ylim(ylim_max * [-1, 1]);
%
% % --- TILE 2: Slope Low Freq (Z-score) ---
% nexttile(t);
% data_plot = zscore(slope_low);
% bar(1:num_ics, zscore(slope_low), 'FaceColor', bar_color, 'EdgeColor', 'none', 'BarWidth', 0.6);
% grid on;
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Box', 'off', 'FontName', font_name, 'FontSize', 10);
% title('Low Frequency Slope (Z-Score)', 'FontSize', 11, 'FontWeight', 'bold');
% xlim([0.5, num_ics + 0.5]);
% xticks(1:num_ics);
% xticklabels({});
% ylim_max = max(4, round(max(abs(data_plot))));
% ylim(ylim_max * [-1, 1]);
%
% % --- TILE 3: Slope Discrepancy (Z-score) ---
% nexttile(t);
% data_plot = zscore(slope_discrepancy);
% bar(1:num_ics, zscore(slope_discrepancy), 'FaceColor', bar_color, 'EdgeColor', 'none', 'BarWidth', 0.6);
% grid on;
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Box', 'off', 'FontName', font_name, 'FontSize', 10);
% title('Slope Discrepancy (Z-Score)', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('Independent Components (ICs)', 'FontWeight', 'bold');
% ylabel('Z-Score');
% xlim([0.5, num_ics + 0.5]);
% xticks(1:num_ics);
% ylim_max = max(4, round(max(abs(data_plot))));
% ylim(ylim_max * [-1, 1]);
%
% % --- TILE 4: Final Normalised Ratio (r) ---
% nexttile(t);
% bar(1:num_ics, r, 'FaceColor', ratio_color, 'EdgeColor', 'none', 'BarWidth', 0.6);
% hold on;
%
% % Add a horizontal dashed line for your target EMG threshold
% line_thresh = yline(threshold_emg_r, '--', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.5);
%
% grid on;
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Box', 'off', 'FontName', font_name, 'FontSize', 10);
% title('Normalised EMG Detection Ratio (r)', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('Independent Components (ICs)', 'FontWeight', 'bold');
% ylabel('Ratio Value');
% xlim([0.5, num_ics + 0.5]);
% ylim_max = max(10, round(max(r)));
% ylim([0, ylim_max]);
% xticks(1:num_ics);
%
% % Add a clean, unobtrusive legend for the decision threshold
% legend(line_thresh, 'EMG Threshold', 'Location', 'NorthEast', 'Box', 'off');
%
% % Link all x-axes across the 2x2 grid for synchronized zooming and panning
% linkaxes(findall(fh, 'type', 'axes'), 'x');
%
% plotX = 30;
% plotY = plotX / 1.6;
% save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_emg'], [plotX plotY]);

%% ========================================================================
% 5. Use icablinkmetrics for eyeblinks
fprintf('\n--------------------------------\n');
fprintf('VEOG: icablinkmetrics plugin\n');
fprintf('--------------------------------\n');

% Blink / VEOG channel
% data_eog_blink = mean(EEG.data(ismember({EXT.chanlocs.labels}, 'VEOG'), :), 1); % fails for some reason
data_eog_blink = mean(EEG.data(ismember({EEG.chanlocs.labels}, cfg.ica.blinkchans), :), 1);

% Filter EEG (~EOG)
% data_eog_blink = do_filteringcore(blEOG,alEOG,EOGdata,EEG.event,EEG.srate);
% data_eog_blink = do_filteringcore(bhEOG,ahEOG,EOGdata,EEG.event,EEG.srate);
% data_eog_blink = data_eog_blink(:, 2)';

% Put bandpassed ICA data
% EEG.icaact = data_ica_eog';
EEG.icaact = data_ica;

try
    icablinkmetricsout = icablinkmetrics(EEG, 'ArtifactChannel', data_eog_blink, 'Alpha', 0.001, 'VisualizeData', 'False');
    if any(icablinkmetricsout.identifiedcomponents > 0)
        fprintf('Blink ICs (N = %d) were identified using icablinkmetrics.\n', length(icablinkmetricsout.identifiedcomponents));
        str = strjoin(string(icablinkmetricsout.identifiedcomponents), ', ');
        fprintf('VEOG ICs: %s\n', str);
    else
        fprintf('No blink ICs were identified using icablinkmetrics.\n');
        icablinkmetricsout.identifiedcomponents = [];
        % icablinkmetricsout.metrics.corr_Pvalue  = [];
        % icablinkmetricsout.metrics.conv_Pvalue  = [];
        % icablinkmetricsout.metrics.perc_Pvalue  = [];
    end
catch
    fprintf('The method has failed. Skipping...\n');
    icablinkmetricsout.identifiedcomponents = [];
    icablinkmetricsout.metrics.corr_Pvalue  = [];
    icablinkmetricsout.metrics.conv_Pvalue  = [];
    icablinkmetricsout.metrics.perc_Pvalue  = [];
end

% Remove
EEG.icaact = [];

% Log
EEG.ALSUTRECHT.ica.blink.icablinkmetrics.pval = single([icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]');
EEG.ALSUTRECHT.ica.blink.icablinkmetrics.bics = icablinkmetricsout.identifiedcomponents(:);
EEG.ALSUTRECHT.ica.blink.icablinkmetrics.cvec = ones(length(EEG.ALSUTRECHT.ica.blink.icablinkmetrics.bics),1);
EEG.ALSUTRECHT.ica.blink.icablinkmetrics.clss = 'Blink';

%% ========================================================================
% 6. Use spatial characteristics for channel pops and wobbles
fprintf('\n--------------------------------\n');
fprintf('Channel: Spatial characteristics\n');
fprintf('--------------------------------\n');

% Estimate spatial smootheness
[spatialSmoothness, bad_ic] = estimate_spatialsmoothnes(EEG, threshold_channel_kurt);

% % These are likely not good/clear ICs (or muscle)
% % -> multiple strong weights
% falseChannel = abs(zscore(icawinv(:,badIC)));
% ICsMostLikelyChannelWrong = sum(falseChannel > 3, 1) ~= 1;
% badIC(ICsMostLikelyChannelWrong) = [];

% falseChannel = sign(falseChannel) .* falseChannel.^2;
% maxWeights = max(abs(falseChannel));

% NICAbad = length(badIC);
% ICsMostLikelyChannelWrong1 = false(1,NICAbad);
% ICsMostLikelyChannelWrong2 = false(1,NICAbad);
%
% for i = 1:NICAbad
%     maxWeightsTmp1 = 0.5 * maxWeights(i);
%     maxWeightsTmp2 = 0.2 * maxWeights(i);
%     ICsMostLikelyChannelWrong1(i) = sum(falseChannel(:,i)>maxWeightsTmp1) > 0 & sum(falseChannel(:,i)<-maxWeightsTmp1) > 0; % EMG dipoles
%     ICsMostLikelyChannelWrong2(i) = sum(abs(falseChannel(:,i))>maxWeightsTmp2) > 3; % too weak channel ICs
% end
% badIC(ICsMostLikelyChannelWrong1 | ICsMostLikelyChannelWrong2) = [];
% figure;
% for i = 1:length(badIC)
%     mytopoplot(icawinv(:,badIC(i)), [],num2str(badIC(i)),nexttile);
% end

% Report
% report_detections('Channel', 'spatial characteristics', {'Cutoff1', 'Cutoff2'}, badIC, spatialTreshold);
report_detections('Channel', 'spatial characteristics', 'kurtosis', bad_ic, threshold_channel_kurt);

% Log
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.pval = spatialSmoothness;
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.bics = bad_ic(:);
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.cvec = ones(length(EEG.ALSUTRECHT.ica.channel.spatialSmoothness.bics),1);
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.clss = 'Channel';

%% ========================================================================
% 7. Use blink IC template for blinks
fprintf('\n--------------------------------\n');
fprintf('VEOG: Blink IC template\n');
fprintf('--------------------------------\n');

% % Lower to capture imperfect ICs, but leads to false positives
% corrMat1 = abs(corr(icawinvSmooth, templates_ica.Blinkweights0, "type", "Spearman"));
% corrMat2 = abs(corr(icawinvSmooth, templates_ica.Blinkweights1, "type", "Spearman"));
%
% badIC1 = find(corrMat1 > threshold_temp_blink);
% badIC2 = find(corrMat2 > threshold_temp_blink);
% bad_ic  = unique([badIC1(:); badIC2(:)]);
%
% Log
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.pval1 = corrMat1;
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.pval2 = corrMat2;
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics1 = badIC1(:);
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics2 = badIC2(:);
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics  = bad_ic(:);
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.cvec  = ones(length(EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics),1);
% EEG.ALSUTRECHT.ica.blink.TemplateCorr.clss  = 'Blink';

% New
cfg_base = struct('check_distribution', false, 'sim_thresh', threshold_temp_blink);
[bad_ic, ~, ~, report] = match_ica_template(icawinv, EEG.chanlocs, templates_ica, 'blink', cfg_base);

% Report
report_detections('VEOG', 'template similarity', 'S', bad_ic, threshold_temp_blink);

% Log
EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics   = bad_ic(:);
EEG.ALSUTRECHT.ica.blink.TemplateCorr.report = report;

%% ========================================================================
% 8. Use saccade IC template
fprintf('\n--------------------------------\n');
fprintf('HEOG: Saccade IC template\n');
fprintf('--------------------------------\n');

% % Lower to capture imperfect ICs, but leads to false positives
% corrMat1 = abs(corr(icawinvSmooth, templates_ica.Saccadeweights0,  "type", "Spearman"));
% corrMat2 = abs(corr(icawinvSmooth, templates_ica.Saccadeweights2L, "type", "Spearman"));
% corrMat3 = abs(corr(icawinvSmooth, templates_ica.Saccadeweights2R, "type", "Spearman"));
%
% badIC1 = find(corrMat1 > threshold_temp_saccade);
% badIC2 = find(corrMat2 > threshold_temp_saccade);
% badIC3 = find(corrMat3 > threshold_temp_saccade);
% bad_ic  = unique([badIC1(:); badIC2(:); badIC3(:)]);
%
% % Report
% report_detections('HEOG', 'template correlation', 'R', badIC1, threshold_temp_saccade);
% report_detections('HEOG', 'template correlation', 'R', badIC2, threshold_temp_saccade);
% report_detections('HEOG', 'template correlation', 'R', badIC3, threshold_temp_saccade);
%
% % Log
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval1 = corrMat1;
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval2 = corrMat2;
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval3 = corrMat3;
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics1 = badIC1(:);
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics2 = badIC2(:);
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics3 = badIC3(:);
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics  = bad_ic(:);
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.cvec  = ones(length(EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics),1);
% EEG.ALSUTRECHT.ica.saccade.TemplateCorr.clss  = 'Saccade';

% New
cfg_base = struct('check_distribution', false, 'sim_thresh', threshold_temp_saccade);
[bad_ic, ~, ~, report] = match_ica_template(icawinv, EEG.chanlocs, templates_ica, 'saccade', cfg_base);

% Report
report_detections('HEOG', 'template similarity', 'S', bad_ic, threshold_temp_saccade);

% Log
EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics   = bad_ic(:);
EEG.ALSUTRECHT.ica.saccade.TemplateCorr.report = report;

%% ========================================================================
% 9. Use heart IC template
fprintf('\n--------------------------------\n');
fprintf('ECG: Heart IC template\n');
fprintf('--------------------------------\n');

% % Lower to capture imperfect ICs, but leads to false positives
% corrMat1 = abs(corr(icawinv, templates_ica.Heartweights0, "type", "Spearman"));
% corrMat2 = abs(corr(icawinv, templates_ica.Heartweights1, "type", "Spearman"));
% corrMat3 = abs(corr(icawinv, templates_ica.Heartweights2, "type", "Spearman"));
%
% badIC1 = find(any(corrMat1 > threshold_temp_heart, 2));
% badIC2 = find(corrMat2 > threshold_temp_heart);
% badIC3 = find(corrMat3 > threshold_temp_heart);
% bad_ic  = unique([badIC1(:); badIC2(:); badIC3(:)]);
%
% % Report
% report_detections('ECG', 'template correlation', 'R', badIC1, threshold_temp_heart);
% report_detections('ECG', 'template correlation', 'R', badIC2, threshold_temp_heart);
% report_detections('ECG', 'template correlation', 'R', badIC3, threshold_temp_heart);
%
% % Log
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval1 = corrMat1;
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval2 = corrMat2;
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval3 = corrMat3;
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics1 = badIC1(:);
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics2 = badIC2(:);
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics3 = badIC3(:);
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics  = bad_ic(:);
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.cvec  = ones(length(EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics), 1);
% EEG.ALSUTRECHT.ica.heart.TemplateCorr.clss  = 'Heart';

% New
cfg_base = struct('check_distribution', false, 'sim_thresh', threshold_temp_heart);
[bad_ic, ~, ~, report] = match_ica_template(icawinv, EEG.chanlocs, templates_ica, 'heart', cfg_base);

% Report
report_detections('ECG', 'template similarity', 'S', bad_ic, threshold_temp_heart);

% Log
EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics   = bad_ic(:);
EEG.ALSUTRECHT.ica.heart.TemplateCorr.report = report;

%% ========================================================================
% 10. Bad ICs
fprintf('\n--------------------------------\n');
fprintf('Bad (general): Spectral power (negative integral)\n');
fprintf('--------------------------------\n');

% Check power spectra properties
[neg_integral, slope_low, slope_discrepancy] = estimate_spectral_metrics(fftBins, freq, [1 70], [2 40]);

badIC1 = find(neg_integral > threshold_power_int);
bad_ic = badIC1(:);

% Report
report_detections('BAD', 'negative integral', 'S', bad_ic, threshold_power_int);

% Log
EEG.ALSUTRECHT.ica.bad.powerspectra.neg_integral = neg_integral;
EEG.ALSUTRECHT.ica.bad.powerspectra.bics1 = badIC1(:);
EEG.ALSUTRECHT.ica.bad.powerspectra.bics  = bad_ic(:);
EEG.ALSUTRECHT.ica.bad.powerspectra.cvec  = zeros(length(EEG.ALSUTRECHT.ica.bad.powerspectra.bics), 1);
EEG.ALSUTRECHT.ica.bad.powerspectra.clss  = 'Bad';

%% ========================================================================
% Final log
fprintf('\n--------------------------------\n');
fprintf('Combining detected ICs\n');
fprintf('--------------------------------\n');

% *Blink ICs
% blink1 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec == 2);
% blink2 = EEG.ALSUTRECHT.ica.blink.icablinkmetrics.bics;
% blink3 = EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics;
% blink = is_false_veog(blink, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);
% ICsMostLikelyBlink = false(num_ics, 1);
% ICsMostLikelyBlink(blink) = true;

m = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec == 2);
blink1 = false(num_ica, 1);
blink1(m) = true;
m = EEG.ALSUTRECHT.ica.blink.icablinkmetrics.bics;
blink2 = false(num_ica, 1);
blink2(m) = true;
m = EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics;
blink3 = false(num_ica, 1);
blink3(m) = true;

blink = [blink1(:), blink2(:), blink3(:)];
blink = mean(blink, 2);
blink = find(blink > 0.5);

blink = is_false_veog(blink, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);
ICsMostLikelyBlink = false(num_ica, 1);
ICsMostLikelyBlink(blink) = true;

% *Saccades ICs
saccade1 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec == 3);
saccade2 = EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics;
saccade  = unique([saccade1(:); saccade2(:)]);
saccade  = is_false_heog(saccade, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);
% if length(saccade) > 1
%     likely_heog = select_likely_heog(saccade, EEG);
%     saccade = saccade(likely_heog);
% end

ICsMostLikelySaccade = false(num_ica, 1);
ICsMostLikelySaccade(saccade) = true;

% *Eye ICs
eye = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 3);
ICsMostLikelyEyeICLabel = false(num_ica, 1);
ICsMostLikelyEyeICLabel(eye) = true;

% % Fix wrong ICLabel estimates
% % 1. Increase the probability, > 0.5/0.6
% % 2. Remove ICs with low variance, likely channel (or muscle) ICs
% ICsMostLikelyEyeICLabel(EEG.ALSUTRECHT.ica.spatialSmoothness.bics) = false;

% Eye ICs: combine all
ICsMostLikelyEye = ICsMostLikelyBlink | ICsMostLikelySaccade | ICsMostLikelyEyeICLabel;

% *Muscle ICs
muscle1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 2);
muscle2 = EEG.ALSUTRECHT.ica.muscle.slope.bics;
muscle  = unique([muscle1(:); muscle2(:)]);
muscle  = is_false_emg(muscle, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);

ICsMostLikelyMuscle = false(num_ica, 1);
ICsMostLikelyMuscle(muscle) = true;

% *Complex ICs
ICsMostLikelyComplex = ICsMostLikelyMuscle & ICsMostLikelyEye;

% Update
ICsMostLikelyMuscle(ICsMostLikelyComplex)     = false;
ICsMostLikelyBlink(ICsMostLikelyComplex)      = false;
ICsMostLikelySaccade(ICsMostLikelyComplex)    = false;
ICsMostLikelyEyeICLabel(ICsMostLikelyComplex) = false;
ICsMostLikelyEye = ICsMostLikelyBlink | ICsMostLikelySaccade | ICsMostLikelyEyeICLabel;

% *Channel ICs
channel1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 6);
channel2 = EEG.ALSUTRECHT.ica.channel.spatialSmoothness.bics;
channel  = unique([channel1(:); channel2(:)]);
ICsMostLikelyChannel = false(num_ica, 1);
ICsMostLikelyChannel(channel) = true;

% Sometimes channel ICs are marked as muscle ICs
% Also ensure that there is no overlap with other IC types
ICsMostLikelyChannelWrong = ICsMostLikelyChannel & (ICsMostLikelyEye | ICsMostLikelyMuscle | ICsMostLikelyComplex);
ICsMostLikelyChannel(ICsMostLikelyChannelWrong) = false;

% *Heart ICs
heart1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 4);
heart2 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec == 1);
heart3 = EEG.ALSUTRECHT.ica.heart.ctps.bics;
heart4 = EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics;
heart  = unique([heart1(:); heart2(:); heart3(:); heart4(:)]);
heart  = is_false_ecg(heart, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);

% % Is the problem correlation with topoplots?
% % I could add here cross-trial phase statistics as well
% if length(heart) > 1
%     fprintf('More than 1 heart IC detected. Removing possible multiple detections based on template correlations.\n');
%     heart_mostlikely = select_likely_ecg(heart4, EEG);
%     heart = unique([heart1(:); heart2(:); heart3(:); heart_mostlikely(:)]);
%     heart = is_false_ecg(heart, icawinv, templates_ica, EEG.ALSUTRECHT.ica.ICLabel);
% end
%
% % If ICLAbel markes it as 'heart' (very rare tho), then it is definitely correct
% heart = unique([heart(:) heart1(:)]);
%
% % Not sure what else to do here
% if length(heart) > 1
%     fprintf('Still more than 1 heart IC detected!\n');
% end

ICsMostLikelyHeart = false(num_ica, 1);
ICsMostLikelyHeart(heart) = true;

if ~isnan(is_ecg)
    ICsMostLikelyHeart = ICsMostLikelyHeart(:) & is_ecg(:);
end

% Ensure no overlap
% ICsMostLikelyHeartWrong = ICsMostLikelyHeart & (ICsMostLikelyEye | ICsMostLikelyMuscle |  ICsMostLikelyComplex | ICsMostLikelyChannel);
% ICsMostLikelyHeart(ICsMostLikelyHeartWrong) = false;
ICsMostLikelyEye(ICsMostLikelyHeart)     = false;
ICsMostLikelyMuscle(ICsMostLikelyHeart)  = false;
ICsMostLikelyComplex(ICsMostLikelyHeart) = false;
ICsMostLikelyChannel(ICsMostLikelyHeart) = false;

% There should be no overlap
assert(max(sum([ICsMostLikelyEye, ICsMostLikelyMuscle, ICsMostLikelyComplex, ICsMostLikelyHeart, ICsMostLikelyChannel], 2)) == 1);

% Bad
bad = EEG.ALSUTRECHT.ica.bad.powerspectra.bics;
false_bad = is_false_brain(bad, EEG.ALSUTRECHT.ica.ICLabel, 0.60);
bad(false_bad) = [];

ICsMostLikelyBad = false(num_ica, 1);
ICsMostLikelyBad(bad) = true;

% ics_bad_sum = any([ICsMostLikelyEye, ICsMostLikelyMuscle, ICsMostLikelyComplex, ICsMostLikelyHeart, ICsMostLikelyChannel], 2);
% x = ICsMostLikelyBad(ics_bad_sum);
% x = x(1:35);

% Final Logging
EEG.ALSUTRECHT.ica.final.eye     = ICsMostLikelyEye;
EEG.ALSUTRECHT.ica.final.muscle  = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.final.complex = ICsMostLikelyComplex;
EEG.ALSUTRECHT.ica.final.channel = ICsMostLikelyChannel;
EEG.ALSUTRECHT.ica.final.heart   = ICsMostLikelyHeart;
EEG.ALSUTRECHT.ica.final.genbad  = ICsMostLikelyBad;

% Labels:
% 1 'Brain'
% 2 'Muscle'
% 3 'Eye'
% 4 'Heart'
% 5 'Line Noise'
% 6 'Channel Noise'
% 7 'Other'
EEG.ALSUTRECHT.ica.final.report = EEG.ALSUTRECHT.ica.ICLabel.cvec;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyMuscle)  = 2;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyEye)     = 3;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyHeart)   = 4;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyChannel) = 6;

% -------------------------------------------------------------------------
% Estimate variance for first K ICs

% Update
EEG = eeg_checkset(EEG, 'ica');
EEG.icaact = [];

% Ensure total number of ICs is defined
num_ica = size(EEG.ALSUTRECHT.ica.icaweights, 1);

% Make sure IC activations are present
ch_idx = EEG.ALSUTRECHT.ica.icachansind;
icaact = reshape(EEG.data(ch_idx, :, :), length(ch_idx), []);
icaact = (EEG.ALSUTRECHT.ica.icaweights * EEG.ALSUTRECHT.ica.icasphere) * icaact;

% Helper function handle for joint VAF calculation
get_class_vaf = @(mask) get_vaf(EEG, icaact, mask);

% The first K ICs
num_ica_relevant = min(20, num_ica);

mask_relevant = false(num_ica, 1);
mask_relevant(1:num_ica_relevant) = true;
var_relevant = get_class_vaf(mask_relevant);

fprintf('\nWithin the first %d ICs (scalp variance accounted for = %.2f%%):\n', num_ica_relevant, var_relevant);

% -------------------------------------------------------------------------
% Class percentage breakdown within first K ICs
cat_labels = {'brain', 'muscle', 'eye', 'heart', 'line', 'channel', 'other'};
report_rel = EEG.ALSUTRECHT.ica.final.report(1:num_ica_relevant);

% Class Variance Accounted For across all ICs
I = EEG.ALSUTRECHT.ica.final.report(:)';
vaf_cats = zeros(1, length(cat_labels));

for i = 1:length(cat_labels)
    pct_comp = mean(report_rel == i) * 100;

    % True scalp VAF (0-100%) for each IC category across all components
    vaf_cats(i) = get_class_vaf(I == i);
    EEG.ALSUTRECHT.ica.final.var.(cat_labels{i}) = vaf_cats(i);

    fprintf('%-10s components: %2.0f%% (scalp variance accounted for = %.2f%%)\n', cat_labels{i}, round(pct_comp), vaf_cats(i));
end

% -------------------------------------------------------------------------
% Handle general bad components explicitly
genbad_mask = EEG.ALSUTRECHT.ica.final.genbad(:)';
vaf_bad = get_class_vaf(genbad_mask);
EEG.ALSUTRECHT.ica.final.var.genbad = vaf_bad;

pct_bad_rel = mean(EEG.ALSUTRECHT.ica.final.genbad(1:num_ica_relevant)) * 100;
fprintf('%-10s components: %2.0f%% (scalp variance accounted for = %.2f%%)\n', 'bad', round(pct_bad_rel), vaf_bad);


end

% function bad_ic = is_false_veog(bad_ic, ICs, templates_ica, ICLabel)
% if isempty(bad_ic); return; end
% false_1 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights0, "type", "Spearman")) < 0.80;
% false_2 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights1, "type", "Spearman")) < 0.80;
% % false_3 = is_false(bad_ic, ICLabel); % can fail tho!
% % assert(isequal(length(bad_ic), length(false_1), length(false_2), length(false_3)));
% assert(isequal(length(bad_ic), length(false_1), length(false_2)));
% false_all = false_1(:) & false_2(:);
% bad_ic(false_all) = [];
% fprintf('Found %d/%d false VEOG ICs. \n', sum(false_all), length(false_all));
% end
%
% function bad_ic = is_false_heog(bad_ic, ICs, templates_ica, ICLabel)
% % VEOG ICs are sometimes very correlated with the HEOG signal
% if isempty(bad_ic); return; end
% false_1 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights0, "type", "Spearman")) > 0.85;
% false_2 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights1, "type", "Spearman")) > 0.85;
% false_3 = is_false(bad_ic, ICLabel);
% assert(isequal(length(bad_ic), length(false_1), length(false_2), length(false_3)));
% false_all = false_1(:) | false_2(:) | false_3(:);
% bad_ic(false_all) = [];
% fprintf('Found %d/%d false HEOG ICs.\n', sum(false_all), length(false_all));
% end
%
% function bad_ic = is_false_emg(bad_ic, ICs, templates_ica, ICLabel)
% if isempty(bad_ic); return; end
% false_1 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights0, "type", "Spearman")) > 0.95;
% false_2 = abs(corr(ICs(:, bad_ic), templates_ica.Blinkweights1, "type", "Spearman")) > 0.95;
% false_3 = is_false(bad_ic, ICLabel);
% assert(isequal(length(bad_ic), length(false_1), length(false_2), length(false_3)));
% false_all = false_1(:) | false_2(:) | false_3(:);
% bad_ic(false_all) = [];
% fprintf('Found %d/%d false EMG ICs.\n', sum(false_all), length(false_all));
% end
%
% function bad_ic = is_false_ecg(bad_ic, ICs, templates_ica, ICLabel)
% if isempty(bad_ic); return; end
% false_1 = abs(corr(ICs(:, bad_ic), templates_ica.Heartweights0, "type", "Pearson")) < 0.75;
% false_1 = all(false_1, 2);
% false_2 = abs(corr(ICs(:, bad_ic), templates_ica.Heartweights1, "type", "Spearman")) < 0.75;
% false_3 = abs(corr(ICs(:, bad_ic), templates_ica.Heartweights2, "type", "Spearman")) < 0.75;
% % % ECG ICs is almsot oalways marked as 'brain' but typically never with a very high confidence
% % false_4 = ICLabel.cvec(bad_ic) == 1 & ICLabel.pvec(bad_ic) >= 0.98;
% false_4 = false(size(bad_ic));
% assert(isequal(length(bad_ic), length(false_1), length(false_2), length(false_3), length(false_4)));
% false_all = (false_1(:) & false_2(:) & false_3(:)) | false_4(:);
% bad_ic(false_all) = [];
% fprintf('Found %d/%d false ECG ICs.\n', sum(false_all), length(false_all));
% end
%
% function false_artifact = is_false(bad_ic, ICLabel)
% % HEOG ICs are sometimes confused by broad L-R dipolar brain ICs
% % EMG ICs are rarely confused by broad almost monopolar brain ICs
% false_artifact = ICLabel.cvec(bad_ic) == 1 & ICLabel.pvec(bad_ic) > 0.6;
% % fprintf('Found %d false artifact ICs.\n', sum(false_all));
% % bad_ic(false_artifact) = [];
% end
% function bad_ic_final = select_likely_heog(bad_ic, EEG)
% % max R -> the most likely IC; There is typically only one
% % bad_ic = EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics;
% r_1 = EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval1(bad_ic, :);
% r_2 = EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval2(bad_ic, :);
% r_3 = EEG.ALSUTRECHT.ica.saccade.TemplateCorr.pval3(bad_ic, :);
% r   = [r_1, r_2, r_3];
% [p, mostlikely] = max(max(r, [], 2));
% bad_ic_final = bad_ic(mostlikely);
% fprintf('From %d HEOG ICs, the most likely HEOG IC is %d.\n', length(bad_ic), bad_ic_final);
% end

% function bad_ic_final = select_likely_ecg(bad_ic, EEG)
% % bad_ic = EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics;
% r_1 = EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval1(bad_ic, :);
% r_2 = EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval2(bad_ic, :);
% r_3 = EEG.ALSUTRECHT.ica.heart.TemplateCorr.pval3(bad_ic, :);
% r   = [r_1, r_2, r_3];
% % max R -> the most likely ECG IC
% % There is typically only one ECG IC
% [p, mostlikely] = max(max(r, [], 2));
% bad_ic_final = bad_ic(mostlikely);
% fprintf('From %d ECG ICs, the most likely ECG IC is %d.\n', length(bad_ic), bad_ic_final);
% end

% =========================================================================
% FALSE-POSITIVE VALIDATION SUITE FOR CANDIDATE ARTIFACT ICs
% =========================================================================

function bad_ic = is_false_veog(bad_ic, ICs, templates_ica, ICLabel)
% Removes candidate VEOG ICs that fail to match blink templates
% or are strongly classified as cortical brain activity.
if isempty(bad_ic); return; end

% Gather all blink template variants (e.g. Blinkweights0 [128x1], Blinkweights1 [128x1])
blink_templates = [templates_ica.Blinkweights0, templates_ica.Blinkweights1];

% Maximum similarity across all blink templates [1 x length(bad_ic)]
max_sim_blink = calc_max_sim(ICs(:, bad_ic), blink_templates, 'cosine');

% Condition 1: Fails to match any blink template (< 0.80)
false_no_template = max_sim_blink < 0.80;

% Condition 2: Confirmed cortical brain source by ICLabel (P > 0.70)
% false_is_brain = is_false_brain(bad_ic, ICLabel, 0.70); % it can also fail (if not a lot of blinks are present)
false_is_brain = false(size(false_no_template));

% Flag as false positive if it fails template matching OR is cortical brain
false_all = false_no_template(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[VEOG Validation] Pruned %d/%d false-positive VEOG ICs.\n', sum(false_all), length(false_all));
end


function bad_ic = is_false_heog(bad_ic, ICs, templates_ica, ICLabel)
% Removes candidate HEOG ICs that are actually blinks (VEOG) or cortical brain.
% VEOG ICA activations are sometimes very correlated with the recorded HEOG signal
if isempty(bad_ic); return; end

% Gather blink templates to catch vertical bleed-through
blink_templates = [templates_ica.Blinkweights0, templates_ica.Blinkweights1];
max_sim_blink = calc_max_sim(ICs(:, bad_ic), blink_templates, 'cosine');

% Condition 1: Component is actually a blink misclassified as HEOG (> 0.85)
false_is_blink = max_sim_blink > 0.85;

% Condition 2: Confirmed cortical brain source (e.g. lateral temporal beta/mu, P > 0.60)
false_is_brain = is_false_brain(bad_ic, ICLabel, 0.60);

false_all = false_is_blink(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[HEOG Validation] Pruned %d/%d false-positive HEOG ICs.\n', sum(false_all), length(false_all));
end


function bad_ic = is_false_emg(bad_ic, ICs, templates_ica, ICLabel)
% =========================================================================
% IS_FALSE_EMG: Prunes false-positive EMG candidate ICs
% =========================================================================
% Removes candidate EMG ICs that are actually:
%   1. Ocular blinks / VEOG (cosine similarity > 0.90)
%   2. Horizontal saccades / HEOG (cosine similarity > 0.85)
%   3. Genuine cortical brain rhythms (ICLabel P(Brain) >= 0.60)
% =========================================================================

if isempty(bad_ic); return; end

% Blink / VEOG Check
blink_fields = fieldnames(templates_ica);
blink_mask   = contains(lower(blink_fields), 'blink') | contains(lower(blink_fields), 'vertical');
blink_templates = [];
for i_f = find(blink_mask)'
    blink_templates = [blink_templates, double(templates_ica.(blink_fields{i_f}))]; %#ok<AGROW>
end

if ~isempty(blink_templates)
    max_sim_blink = calc_max_sim(ICs(:, bad_ic), blink_templates, 'cosine');
    false_is_blink = max_sim_blink > 0.90;
else
    false_is_blink = false(1, length(bad_ic));
end

% Saccade / HEOG Check
heog_mask = contains(lower(blink_fields), 'horizontal') | ...
    contains(lower(blink_fields), 'saccade')    | ...
    contains(lower(blink_fields), 'heog');
heog_templates = [];
for i_f = find(heog_mask)'
    heog_templates = [heog_templates, double(templates_ica.(blink_fields{i_f}))]; %#ok<AGROW>
end

if ~isempty(heog_templates)
    max_sim_heog = calc_max_sim(ICs(:, bad_ic), heog_templates, 'cosine');
    false_is_heog = max_sim_heog > 0.85;
else
    false_is_heog = false(1, length(bad_ic));
end

% Cortical Brain Protection Check (ICLabel)
% Protects genuine sensorimotor beta/gamma or focal temporal alpha
false_is_brain = is_false_brain(bad_ic, ICLabel, 0.60);

% Combine Rejections
false_all = false_is_blink(:) | false_is_heog(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[EMG Validation]  Pruned %d/%d false-positive EMG ICs (Blink: %d, HEOG: %d, Brain: %d).\n', ...
    sum(false_all), length(false_all), sum(false_is_blink), sum(false_is_heog), sum(false_is_brain));

end


function bad_ic = is_false_ecg(bad_ic, ICs, templates_ica, ICLabel)
% Removes candidate ECG ICs that fail to match any cardiac template variant.
if isempty(bad_ic); return; end

% Concatenate all cardiac template variants:
% Heartweights0 [128x3] + Heartweights1 [128x1] + Heartweights2 [128x1] = [128x5]
heart_templates = [double(templates_ica.Heartweights0), ...
    double(templates_ica.Heartweights1), ...
    double(templates_ica.Heartweights2)];

% Maximum similarity across all 5 cardiac variants
max_sim_heart = calc_max_sim(ICs(:, bad_ic), heart_templates, 'cosine');

% Condition 1: Fails to match any cardiac template (< 0.75)
false_no_template = max_sim_heart < 0.75;

% Condition 2: ICLabel is extremely confident it is Brain (P >= 0.95)
% (Cardiac topographies often get weak brain labels (0.5-0.7), so threshold must be strict)
false_is_brain = is_false_brain(bad_ic, ICLabel, 0.98);

false_all = false_no_template(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[ECG Validation]  Pruned %d/%d false-positive ECG ICs.\n', sum(false_all), length(false_all));
end

function is_brain = is_false_brain(bad_ic, ICLabel, prob_threshold)
% Checks whether candidate ICs are classified as Brain by ICLabel
% ICLabel class 1 = 'Brain'
if isempty(ICLabel) || ~isfield(ICLabel, 'cvec') || ~isfield(ICLabel, 'pvec')
    is_brain = false(size(bad_ic));
    return;
end

is_brain = (ICLabel.cvec(bad_ic) == 1) & (ICLabel.pvec(bad_ic) >= prob_threshold);
end


function [bad_ics, sim_scores, dist_scores, report] = match_ica_template(icawinv, chanlocs, templates_ica, template_type, cfg)
% =========================================================================
% MATCH_ICA_TEMPLATE: Multi-variant spatial template matching (BioSemi 128)
% =========================================================================
% Inputs:
%   icawinv       : [128 x num_ics] ICA mixing matrix (EEG.icawinv).
%   chanlocs      : EEGLAB chanlocs struct array (EEG.chanlocs) or [] if unused.
%   templates_ica : Struct containing template weights (e.g. Heartweights0,
%                   Heartweights1, Blinkweights0, Blinkweights1, etc.).
%   template_type : 'heart', 'blink', 'saccade' (or 'horizontal').
%   cfg           : (Optional) Configuration struct:
%                   - cfg.sim_thresh         : Minimum similarity (default: 0.80)
%                   - cfg.metric             : 'cosine' (default) or 'pearson'
%                   - cfg.check_distribution : true / false (default: false)
%                   - cfg.dist_thresh        : Custom distribution cutoff
%
% Outputs:
%   bad_ics     : Row vector of flagged IC indices (e.g. [1, 4])
%   sim_scores  : [1 x num_ics] Best similarity score per IC across variants
%   dist_scores : [1 x num_ics] Focality/diffuseness scores (NaN if disabled)
%   report      : Struct with detailed per-variant matches and diagnostics
% =========================================================================

% -------------------------------------------------------------------------
% 1. Input Parsing and Validation
% -------------------------------------------------------------------------
if nargin < 5; cfg = struct(); end

assert(~isempty(icawinv) && ismatrix(icawinv), 'icawinv must be a non-empty 2D matrix.');
icawinv = double(icawinv);
[num_chans, num_ics] = size(icawinv);

% Parse channel locations for spatial distribution checks
if ~isempty(chanlocs) && isstruct(chanlocs)
    chan_labels = upper({chanlocs.labels});
    has_coords  = isfield(chanlocs, 'X') && ~isempty(chanlocs(1).X);
else
    chan_labels = {};
    has_coords  = false;
end

% -------------------------------------------------------------------------
% 2. Dynamic Template Extraction from templates_ica
% -------------------------------------------------------------------------
all_field_names = fieldnames(templates_ica);

% Match fields containing the template_type keyword (case-insensitive)
matched_idx = contains(lower(all_field_names), lower(template_type));

% Fallback for saccades if labelled as horizontal
if ~any(matched_idx) && strcmpi(template_type, 'saccade')
    matched_idx = contains(lower(all_field_names), 'horizontal');
end

assert(any(matched_idx), ...
    'No fields matching "%s" found in templates_ica struct.', template_type);

matched_fields = all_field_names(matched_idx);

% Horizontally concatenate all matching template columns into [128 x N_variants]
W_tgt = [];
for i_f = 1:length(matched_fields)
    current_template = double(templates_ica.(matched_fields{i_f}));
    assert(size(current_template, 1) == num_chans, ...
        'Field %s has %d channels; expected %d.', ...
        matched_fields{i_f}, size(current_template, 1), num_chans);
    W_tgt = [W_tgt, current_template]; %#ok<AGROW>
end

num_variants = size(W_tgt, 2);

% -------------------------------------------------------------------------
% 3. Parameter Defaults
% -------------------------------------------------------------------------
if ~isfield(cfg, 'metric');             cfg.metric             = 'cosine'; end
if ~isfield(cfg, 'check_distribution'); cfg.check_distribution = false;    end

if ~isfield(cfg, 'sim_thresh') || isempty(cfg.sim_thresh)
    switch lower(template_type)
        case 'blink',   cfg.sim_thresh = 0.80;
        case 'saccade', cfg.sim_thresh = 0.78;
        case 'heart',   cfg.sim_thresh = 0.75;
        otherwise,      cfg.sim_thresh = 0.80;
    end
end

% -------------------------------------------------------------------------
% 4. Multi-Variant Spatial Similarity Calculation via calc_max_sim
% -------------------------------------------------------------------------
[sim_scores, sim_matrix, best_variant_idx] = calc_max_sim(icawinv, W_tgt, cfg.metric);

% -------------------------------------------------------------------------
% 5. Optional BioSemi 128 Spatial Distribution Check
% -------------------------------------------------------------------------
dist_scores = NaN(1, num_ics);
pass_dist   = true(1, num_ics);

if cfg.check_distribution
    total_energy = sum(icawinv.^2, 1) + eps;

    switch lower(template_type)
        case 'blink'
            if ~isfield(cfg, 'dist_thresh') || isempty(cfg.dist_thresh)
                cfg.dist_thresh = 0.60; % >= 60% power in anterior pole
            end

            % BioSemi 128 Prefrontal / Polar ROI (C-bank anterior perimeter)
            biosemi_blink = {'C15', 'C16', 'C17', 'C18', 'C19', 'C20', ...
                'C21', 'C22', 'C25', 'C26', 'C27', 'C28', ...
                'C29', 'C30', 'C31', 'C32'};
            roi_mask = ismember(chan_labels, biosemi_blink);

            if sum(roi_mask) < 4 && has_coords
                y_coords = [chanlocs.Y];
                roi_mask = y_coords > (0.65 * max(y_coords));
            end

            if any(roi_mask)
                dist_scores = sum(icawinv(roi_mask, :).^2, 1) ./ total_energy;
                pass_dist   = dist_scores >= cfg.dist_thresh;
            else
                warning('Prefrontal BioSemi channels not found; bypassing blink distribution check.');
            end

        case {'saccade', 'horizontal'}
            if ~isfield(cfg, 'dist_thresh') || isempty(cfg.dist_thresh)
                cfg.dist_thresh = 0.50; % >= 50% power in outer lateral leads
            end

            % BioSemi 128 Outer Lateral Temporal/Orbital leads
            biosemi_saccade = {'B19', 'B20', 'B21', 'B22', 'B23', 'B24', 'B25', 'B26', ...
                'D11', 'D12', 'D13', 'D14', 'D19', 'D20', 'D21', 'D22', ...
                'D23', 'D24', 'D25', 'D26', 'C11', 'C12', 'C13', 'C14'};
            roi_mask = ismember(chan_labels, biosemi_saccade);

            if sum(roi_mask) < 4 && has_coords
                x_coords = [chanlocs.X];
                y_coords = [chanlocs.Y];
                roi_mask = abs(x_coords) > (0.60 * max(abs(x_coords))) & (y_coords > -0.2 * max(abs(y_coords)));
            end

            if any(roi_mask)
                dist_scores = sum(icawinv(roi_mask, :).^2, 1) ./ total_energy;
                pass_dist   = dist_scores >= cfg.dist_thresh;
            else
                warning('Lateral BioSemi channels not found; bypassing saccade distribution check.');
            end

        case 'heart'
            if ~isfield(cfg, 'dist_thresh') || isempty(cfg.dist_thresh)
                cfg.dist_thresh = 0.15; % Max single lead must not exceed 15% of total cap power
            end

            % True cardiac is broad/diffuse; rejects sharp focal dipoles
            dist_scores = max(icawinv.^2, [], 1) ./ total_energy;
            pass_dist   = dist_scores <= cfg.dist_thresh;
    end
else
    cfg.dist_thresh = [];
end

% -------------------------------------------------------------------------
% 6. Output Packaging & Logging
% -------------------------------------------------------------------------
pass_sim = sim_scores >= cfg.sim_thresh;
bad_mask = pass_sim & pass_dist;
bad_ics  = find(bad_mask);

report = struct();
report.template_type       = template_type;
report.extracted_fields    = matched_fields;
report.num_variants        = num_variants;
report.metric              = cfg.metric;
report.sim_threshold       = cfg.sim_thresh;
report.check_distribution  = cfg.check_distribution;
report.dist_threshold      = cfg.dist_thresh;
report.sim_matrix          = sim_matrix;
report.sim_scores          = sim_scores;
report.best_variant_idx    = best_variant_idx;
report.dist_scores         = dist_scores;
report.bad_ics             = bad_ics;
report.num_detected        = length(bad_ics);

if cfg.check_distribution
    fprintf('[BioSemi-128 Template: %-7s] Detected %d ICs (Sim >= %.2f | %d variants tested | Dist Check Active)\n', ...
        upper(template_type), length(bad_ics), cfg.sim_thresh, num_variants);
else
    fprintf('[BioSemi-128 Template: %-7s] Detected %d ICs (Sim >= %.2f | %d variants tested | Pure Spatial Match)\n', ...
        upper(template_type), length(bad_ics), cfg.sim_thresh, num_variants);
end
end

% =========================================================================
% HELPER FUNCTION: calc_max_sim
% =========================================================================
function [max_sim, sim_matrix, best_variant_idx] = calc_max_sim(W_ic, W_tgt, metric)
% Computes similarity across all candidate ICs and template variants.
%
% Inputs:
%   W_ic   : [num_chans x num_ics] Candidate IC mixing weights
%   W_tgt  : [num_chans x num_variants] Spatial template columns
%   metric : 'cosine' (default) or 'pearson'
%
% Outputs:
%   max_sim          : [1 x num_ics] Highest similarity score per IC
%   sim_matrix       : [num_variants x num_ics] Full pairwise similarity matrix
%   best_variant_idx : [1 x num_ics] Index of the best matching template variant

if nargin < 3 || isempty(metric); metric = 'cosine'; end

W_ic  = double(W_ic);
W_tgt = double(W_tgt);

switch lower(metric)
    case 'cosine'
        % Uncentred Absolute Cosine Similarity: |u' * v| / (||u|| * ||v||)
        norm_ic   = sqrt(sum(W_ic.^2, 1));               % [1 x num_ics]
        norm_tgt  = sqrt(sum(W_tgt.^2, 1))';             % [num_variants x 1]
        dot_prods = abs(W_tgt' * W_ic);                  % [num_variants x num_ics]

        sim_matrix = dot_prods ./ (norm_tgt * norm_ic + eps);

    case 'pearson'
        % Centred Pearson correlation magnitude
        sim_matrix = abs(corr(W_ic, W_tgt, 'type', 'Pearson'))'; % [num_variants x num_ics]

    otherwise
        error('Unknown metric "%s". Use ''cosine'' or ''pearson''.', metric);
end

[max_sim, best_variant_idx] = max(sim_matrix, [], 1);
end

function report_detections(typeIC, typeMethod, badMetric, bad_ic, badTreshold)
if isscalar(badTreshold)
    badCutoff = sprintf('(%s > %1.2f)', badMetric, badTreshold);
else
    badCutoff = sprintf('(%s < %1.1f & %s < %1.1f)', badMetric{1}, badTreshold(1), badMetric{2}, badTreshold(1));
end
if ~isempty(bad_ic)
    fprintf('%s ICs (N = %d) were identified using %s %s.\n', typeIC, length(bad_ic), typeMethod, badCutoff);
    badICstr = strjoin(string(bad_ic),', ');
    fprintf('%s ICs: %s\n', typeIC, badICstr);
else
    fprintf('No %s ICs were identified using %s %s.\n', typeIC, typeMethod, badCutoff);
end
end

function vaf = get_vaf(EEG, icaact, mask)
idx = find(mask);
if isempty(idx)
    vaf = 0;
else
    ch_data = reshape(EEG.data(EEG.ALSUTRECHT.ica.icachansind, :, :), length(EEG.ALSUTRECHT.ica.icachansind), []);
    [~, vaf] = compvar(ch_data, icaact, EEG.ALSUTRECHT.ica.icawinv, idx);
end
end
