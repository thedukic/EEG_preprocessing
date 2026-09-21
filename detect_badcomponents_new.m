function EEG = detect_badcomponents_new(EEG, EXT, EMG, cfg)
% DETECT_BADCOMPONENTS Detects artifactual ICs by gathering independent
% evidence from multiple metrics, visualising diagnostic criteria, and
% synthesising classifications into a unified audit dashboard.

% 1. Initialisation
fprintf('\n================================\n');
fprintf('Detecting bad ICs\n');
fprintf('================================\n');

% Ensure ICA activations
EEG = eeg_checkset(EEG, 'ica');
EEG.icaact = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :);

% Extract ICA templates
templates_ica = load_ictemplateweights(EEG);

% Initialise evidence structure
evidence = struct();

% 2. Gather Evidence
[EEG, evidence] = gather_iclabel(EEG, cfg, evidence);
[EEG, evidence] = gather_ext_correlations(EEG, EXT, cfg, evidence);
[EEG, evidence] = gather_ecg_ctps(EEG, EXT, cfg, evidence);
[EEG, evidence] = gather_ecg_erp(EEG, EXT, cfg, evidence);
[EEG, evidence] = gather_emg_slopes(EEG, cfg, evidence);
[EEG, evidence] = gather_blinkmetrics(EEG, cfg, evidence);
[EEG, evidence] = gather_spatial_characteristics(EEG, cfg, evidence);
[EEG, evidence] = gather_template_matches(EEG, templates_ica, evidence);
[EEG, evidence] = gather_spectral_metrics(EEG, cfg, evidence);
% [EEG, evidence] = gather_temporal_metrics(EEG, cfg, evidence);

% 3. Synthesise Evidence
EEG = synthesise_evidence(EEG, evidence, templates_ica, cfg);

% 4. Variance Accounted For (VAF) Report and Unified Dashboard
EEG = compute_vaf_report(EEG, evidence, cfg);

% Log
EEG.ALSUTRECHT.ica.evidence = evidence;

end

%% ========================================================================
% EVIDENCE GATHERING HELPER FUNCTIONS
% ========================================================================

function [EEG, evidence] = gather_iclabel(EEG, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('ICLabel\n');
fprintf('--------------------------------\n');

active_thresholds = cfg.ica.iclabel.base;
switch lower(EEG.ALSUTRECHT.subject.task)
    case {'mmn', 'sart'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.erp;
    case {'mt'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.mt;
    case {'rs', 'eo', 'ec'}
        active_thresholds(2, 1) = cfg.ica.iclabel.muscle_thresh.rs;
end

EEG = iclabel(EEG);
EEG = pop_icflag(EEG, active_thresholds);

EEG.ALSUTRECHT.ica.ICLabel.bics    = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel.classes = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);

% Store evidence masks
% IGNORES PROBABILITIES !!!
evidence.iclabel_brain   = (EEG.ALSUTRECHT.ica.ICLabel.cvec == 1);
evidence.iclabel_muscle  = (EEG.ALSUTRECHT.ica.ICLabel.cvec == 2);
evidence.iclabel_eye     = (EEG.ALSUTRECHT.ica.ICLabel.cvec == 3);
evidence.iclabel_heart   = (EEG.ALSUTRECHT.ica.ICLabel.cvec == 4);
evidence.iclabel_channel = (EEG.ALSUTRECHT.ica.ICLabel.cvec == 6);
end

function [EEG, evidence] = gather_ext_correlations(EEG, EXT, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('ECG/VEOG/HEOG: Correlations with the external electrodes\n');
fprintf('--------------------------------\n');

threshold_ext = 3;

channel_ecg  = find(strcmp({EXT.chanlocs.labels}, 'ECG'));
channel_veog = find(strcmp({EXT.chanlocs.labels}, 'VEOG'));
channel_heog = find(strcmp({EXT.chanlocs.labels}, 'HEOG'));

% Initialise masks
data_ica = EEG.icaact;
num_ica = size(data_ica, 1);
evidence.corr_veog = false(num_ica, 1);
evidence.corr_heog = false(num_ica, 1);
evidence.corr_ecg  = false(num_ica, 1);

if isempty(channel_ecg)
    fprintf('Warning: ECG signal was not recorded.\n');
    data_ext_ecg = [];
else
    data_ext_ecg = EXT.data(channel_ecg, :);
end

data_ext_eog = EXT.data([channel_veog, channel_heog], :);

% Filter configurations
[bh_eog, ah_eog] = butter(2, 0.3/(EEG.srate/2), 'high');
[bl_eog, al_eog] = butter(2, 10/(EEG.srate/2),  'low');

data_ica_eog = do_filteringcore(bh_eog, ah_eog, do_filteringcore(bl_eog, al_eog, data_ica, EEG.event, EEG.srate), EEG.event, EEG.srate)';
data_ext_eog = do_filteringcore(bh_eog, ah_eog, do_filteringcore(bl_eog, al_eog, data_ext_eog, EEG.event, EEG.srate), EEG.event, EEG.srate);

if ~isempty(data_ext_ecg)
    [bh_ecg, ah_ecg] = butter(2, 10/(EEG.srate/2), 'high');
    [bl_ecg, al_ecg] = butter(2, 20/(EEG.srate/2), 'low');
    data_ica_ecg = do_filteringcore(bh_ecg, ah_ecg, do_filteringcore(bl_ecg, al_ecg, data_ica, EEG.event, EEG.srate), EEG.event, EEG.srate)';
    data_ext_ecg = do_filteringcore(bh_ecg, ah_ecg, do_filteringcore(bl_ecg, al_ecg, data_ext_ecg, EEG.event, EEG.srate), EEG.event, EEG.srate);
    corr_ecg     = corr(data_ica_ecg.^2, abs(data_ext_ecg'), "type", "Spearman");
else
    corr_ecg = NaN(num_ica, 1);
end

corr_veog = corr(data_ica_eog, data_ext_eog(1, :)', "type", "Spearman");
corr_heog = corr(data_ica_eog, data_ext_eog(2, :)', "type", "Spearman");

% Absolute Z-scores across components
corr_ext = abs(zscore([corr_ecg, corr_veog, corr_heog]));

% ECG evidence
if ~any(isnan(corr_ext(:, 1)))
    tmp = find(corr_ext(:, 1) > threshold_ext);
    if length(tmp) > 1
        [~, idx] = max(corr_ext(tmp, 1));
        tmp = tmp(idx);
    end
    evidence.corr_ecg(tmp) = true;
end

% VEOG evidence
evidence.corr_veog(corr_ext(:, 2) > threshold_ext) = true;

% HEOG evidence
tmp = find(corr_ext(:, 3) > threshold_ext);
if length(tmp) > 1
    [~, idx] = max(corr_ext(tmp, 3));
    tmp = tmp(idx);
end
evidence.corr_heog(tmp) = true;

% Log to EEG struct
bad_ic = [find(evidence.corr_ecg); find(evidence.corr_veog); find(evidence.corr_heog)];
bad_ic_type = [ones(sum(evidence.corr_ecg), 1); 2*ones(sum(evidence.corr_veog), 1); 3*ones(sum(evidence.corr_heog), 1)];
EEG.ALSUTRECHT.ica.corr.corr = single(corr_ext);
EEG.ALSUTRECHT.ica.corr.bics = bad_ic;
EEG.ALSUTRECHT.ica.corr.cvec = bad_ic_type;
EEG.ALSUTRECHT.ica.corr.classes = {'ECG', 'VEOG', 'HEOG'};

% Plotting: External Channel Correlations Heatmap
if cfg.figure.plot
    fh = figure('Color', 'w', 'Position', [100, 100, 900, 360], 'Visible', cfg.figure.visible);
    t = tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    ax1 = nexttile(t);

    imagesc(1:num_ica, 1:3, corr_ext');
    colormap(ax1, brewermap([], 'Reds'));
    c = colorbar(ax1);
    c.Label.String = 'Correlation Strength (|Z-Score|)';

    hold on;
    % Mark cells exceeding threshold with a distinct dot
    [flag_ic, flag_ch] = find(corr_ext > threshold_ext);
    if ~isempty(flag_ic)
        plot(ax1, flag_ic, flag_ch, 'w*', 'MarkerSize', 8, 'LineWidth', 1.2);
    end
    hold off;

    set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'YDir', 'normal');
    xlabel('Independent Components (ICs)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Channel', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('External Correlation Heatmap (* indicates |Z| > %g)', threshold_ext), 'FontSize', 12, 'FontWeight', 'bold');
    xticks(1:num_ica);
    yticks(1:3);
    yticklabels({'ECG', 'VEOG', 'HEOG'});

    plotX = max(20, num_ica * 0.7 + 5);
    plotY = 11;
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_external'], [plotX, plotY]);
end
end

function [EEG, evidence] = gather_ecg_ctps(EEG, EXT, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('ECG: Cross-trial phase statistics\n');
fprintf('--------------------------------\n');

threshold_ctps_pk = 20;

data_ica          = EEG.icaact;
num_ica           = size(data_ica, 1);
evidence.ecg_ctps = false(num_ica, 1);
pulse_estimate    = NaN;
V                 = NaN(num_ica, 1);
pK                = NaN(num_ica, 1);
bad_ic            = [];

if any(strcmp({EXT.chanlocs.labels}, 'ECG'))
    cfg_tmp.win          = [-200 250];
    cfg_tmp.do_plot      = true;
    cfg_tmp.plot_visible = cfg.figure.visible;
    [ecg_mask, ecg_epoch, ~, ~, ~, pulse_estimate] = detect_ecg(EXT, cfg_tmp);

    if ~isnan(ecg_mask)
        [V, pK] = my_ctps(data_ica, ecg_epoch, EEG.event, EEG.srate);
        bad_ic = find(pK >= threshold_ctps_pk);
        evidence.ecg_ctps(bad_ic) = true;
    end
end

EEG.ALSUTRECHT.ica.heart.ctps.v    = V;
EEG.ALSUTRECHT.ica.heart.ctps.pk   = pK;
EEG.ALSUTRECHT.ica.heart.ctps.bics = bad_ic(:);
EEG.ALSUTRECHT.subject.heartrate   = pulse_estimate;

% Plotting: CTPS Diagnostics
if cfg.figure.plot && any(~isnan(pK))
    fh = figure('Color', 'w', 'Position', [100, 100, 900, 520], 'Visible', cfg.figure.visible);
    t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

    % Top: pK Significance
    ax1 = nexttile(t);
    b1 = bar(ax1, 1:num_ica, pK, 'FaceColor', [0.35, 0.45, 0.55], 'EdgeColor', 'none', 'BarWidth', 0.6);
    hold(ax1, 'on');
    if any(evidence.ecg_ctps)
        ecg_ctps_tmp = find(evidence.ecg_ctps);
        pK_tmp =  pK(evidence.ecg_ctps);
        for i_comp = 1:length(ecg_ctps_tmp)
            bar(ax1, ecg_ctps_tmp(i_comp), pK_tmp(i_comp), 'FaceColor', [0.85, 0.20, 0.20], 'EdgeColor', 'none', 'BarWidth', b1.BarWidth);
        end
    end
    yline(ax1, threshold_ctps_pk, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.3, 'DisplayName', 'Threshold (p_K \geq 20)');
    hold(ax1, 'off');
    grid(ax1, 'on');
    set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
    ylabel(ax1, 'Significance (p_K)', 'FontSize', 11, 'FontWeight', 'bold');
    title(ax1, 'ICA Phase Locking to ECG R-Peaks (CTPS)', 'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax1, [0.5, num_ica + 0.5]);
    xticks(ax1, 1:num_ica);

    % Bottom: Kuiper Index V
    ax2 = nexttile(t);
    bar(ax2, 1:num_ica, V, 'FaceColor', [0.55, 0.55, 0.55], 'EdgeColor', 'none', 'BarWidth', 0.6);
    grid(ax2, 'on');
    set(ax2, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
    xlabel(ax2, 'Independent Components (ICs)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax2, 'Kuiper Index (V)', 'FontSize', 11, 'FontWeight', 'bold');
    xlim(ax2, [0.5, num_ica + 0.5]);
    ylim(ax2, [0, 1]);
    xticks(ax2, 1:num_ica);

    linkaxes([ax1, ax2], 'x');
    plotX = max(20, num_ica * 0.7 + 5);
    plotY = 14;
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_ctps'], [plotX, plotY]);
end
end

function [EEG, evidence] = gather_ecg_erp(EEG, EXT, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('ECG: Event-related potential alignment\n');
fprintf('--------------------------------\n');

min_corr   = 0.8;
min_snr    = 3.5;
max_lag_ms = 30;

data_ica         = EEG.icaact;
num_ica          = size(data_ica, 1);
evidence.ecg_erp = false(num_ica, 1);
stats            = NaN;
bad_ic           = [];

if any(strcmp({EXT.chanlocs.labels}, 'ECG'))
    cfg_tmp.win          = [-200 250];
    cfg_tmp.do_plot      = true;
    cfg_tmp.plot_visible = cfg.figure.visible;
    [ecg_mask, ecg_epoch] = detect_ecg(EXT, cfg_tmp);

    if ~isnan(ecg_mask)
        cfg_tmp = struct('ecg_chan', 'ECG', 'epoch_window_ms', cfg_tmp.win, 'min_corr', min_corr, 'min_snr', min_snr, ...
            'max_lag_ms', max_lag_ms, 'do_plot', cfg.figure.plot, 'plot_visible', cfg.figure.visible);
        [is_ecg, stats, fh] = check_ic_ecg_erp(EEG, EXT, ecg_epoch, cfg_tmp);

        if cfg.figure.plot && ishandle(fh)
            save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_ecg_erp'], [20, 15]);
        end
        evidence.ecg_erp = is_ecg(:);
        bad_ic = find(is_ecg);
    end
end

EEG.ALSUTRECHT.ica.heart.erp.stats = stats;
EEG.ALSUTRECHT.ica.heart.erp.bics  = bad_ic(:);
end

function [EEG, evidence] = gather_emg_slopes(EEG, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('EMG: Power slopes\n');
fprintf('--------------------------------\n');

data_ica = EEG.icaact;
num_ica  = size(data_ica, 1);
evidence.emg_slope = false(num_ica, 1);

options.Freq_to_compute       = [1 100];
options.muscleFreqEx          = 50 + 2*[-1 1];
options.muscleFreq1           = cfg.emg.slope_freq_1;
options.muscleFreq2           = cfg.emg.slope_freq_2;
options.muscleSlopeThreshold1 = cfg.emg.slope_threshold_1;
options.muscleSlopeThreshold2 = cfg.emg.slope_threshold_2;

[pow, frefAll] = pwelch(data_ica', size(data_ica, 2), [], size(data_ica, 2), EEG.srate);
pow = pow'; frefAll = frefAll';

freq = options.Freq_to_compute(1):0.5:options.Freq_to_compute(2);
fftBins = zeros(size(pow, 1), size(freq, 2));
for i_freq = 1:length(freq)
    [~, index1] = min(abs(frefAll - (freq(i_freq) - 0.25)));
    [~, index2] = min(abs(frefAll - (freq(i_freq) + 0.25)));
    fftBins(:, i_freq) = mean(pow(:, index1:index2), 2);
end

slope_muscle = NaN(num_ica, 2);
for i_ic = 1:num_ica
    % Broad Window
    [~, fin1] = min(abs(options.muscleFreq1(1) - freq));
    [~, fin2] = min(abs(options.muscleFreq1(2) - freq));
    freqHz_broad = freq(fin1:fin2);
    freqPow_broad = fftBins(i_ic, fin1:fin2);

    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_broad));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_broad));
        if fex1 <= length(freqHz_broad)
            freqHz_broad(fex1:fex2) = []; freqPow_broad(fex1:fex2) = [];
        end
    end
    p_broad = polyfit(log10(freqHz_broad), log10(freqPow_broad), 1);

    % High Window
    [~, fhi1] = min(abs(options.muscleFreq2(1) - freq));
    [~, fhi2] = min(abs(options.muscleFreq2(2) - freq));
    freqHz_high = freq(fhi1:fhi2);
    freqPow_high = fftBins(i_ic, fhi1:fhi2);

    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_high));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_high));
        if fex1 <= length(freqHz_high)
            freqHz_high(fex1:fex2) = []; freqPow_high(fex1:fex2) = [];
        end
    end
    p_high = polyfit(log10(freqHz_high), log10(freqPow_high), 1);
    slope_muscle(i_ic, :) = [p_broad(1), p_high(1)];
end

bad_ic = find(slope_muscle(:, 1) > options.muscleSlopeThreshold1 | slope_muscle(:, 2) > options.muscleSlopeThreshold2);
evidence.emg_slope(bad_ic) = true;

EEG.ALSUTRECHT.ica.muscle.slope.slope = slope_muscle;
EEG.ALSUTRECHT.ica.muscle.slope.bics  = bad_ic(:);

% Plotting: Dual-Window EMG Slope Verification
if cfg.figure.plot
    fh = figure('Color', 'w', 'Position', [100, 100, 900, 520], 'Visible', cfg.figure.visible);
    t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

    % Broad Window Slope
    ax1 = nexttile(t);
    bh1 = bar(ax1, 1:num_ica, slope_muscle(:, 1), 'FaceColor', [0.45, 0.45, 0.45], 'EdgeColor', 'none', 'BarWidth', 0.6);
    hold(ax1, 'on');
    flagged_w1 = find(slope_muscle(:, 1) > options.muscleSlopeThreshold1);
    if ~isempty(flagged_w1)
        for i = 1:length(flagged_w1)
            bar(ax1, flagged_w1(i), slope_muscle(flagged_w1(i), 1), 'FaceColor', [0.85, 0.20, 0.20], 'EdgeColor', 'none', 'BarWidth', bh1.BarWidth);
        end
    end
    yline(ax1, options.muscleSlopeThreshold1, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.3);
    hold(ax1, 'off');
    grid(ax1, 'on');
    set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
    ylabel(ax1, 'Slope (Broad)', 'FontSize', 11, 'FontWeight', 'bold');
    title(ax1, sprintf('EMG Dual-Band Spectral Slopes (Broad: %d-%d Hz, High: %d-%d Hz)', ...
        options.muscleFreq1(1), options.muscleFreq1(2), options.muscleFreq2(1), options.muscleFreq2(2)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax1, [0.5, num_ica + 0.5]);
    xticks(ax1, 1:num_ica);

    % High Window Slope
    ax2 = nexttile(t);
    bh2 = bar(ax2, 1:num_ica, slope_muscle(:, 2), 'FaceColor', [0.45, 0.45, 0.45], 'EdgeColor', 'none', 'BarWidth', 0.6);
    hold(ax2, 'on');
    flagged_w2 = find(slope_muscle(:, 2) > options.muscleSlopeThreshold2);
    if ~isempty(flagged_w2)
        for i = 1:length(flagged_w2)
            bar(ax2, flagged_w2(i), slope_muscle(flagged_w2(i), 2), 'FaceColor', [0.85, 0.20, 0.20], 'EdgeColor', 'none', 'BarWidth', bh2.BarWidth);
        end
    end
    yline(ax2, options.muscleSlopeThreshold2, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.3);
    hold(ax2, 'off');
    grid(ax2, 'on');
    set(ax2, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
    xlabel(ax2, 'Independent Components (ICs)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax2, 'Slope (High)', 'FontSize', 11, 'FontWeight', 'bold');
    xlim(ax2, [0.5, num_ica + 0.5]);
    xticks(ax2, 1:num_ica);

    linkaxes([ax1, ax2], 'x');
    plotX = max(20, num_ica * 0.7 + 5);
    plotY = 14;
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_emg_slopes'], [plotX, plotY]);
end
end

function [EEG, evidence] = gather_blinkmetrics(EEG, cfg, evidence)
fprintf('\n--------------------------------\n');
fprintf('VEOG: icablinkmetrics plugin\n');
fprintf('--------------------------------\n');

data_ica              = EEG.icaact;
num_ica               = size(data_ica, 1);
evidence.blinkmetrics = false(num_ica, 1);
data_eog_blink        = mean(EEG.data(ismember({EEG.chanlocs.labels}, cfg.ica.blinkchans), :), 1);

try
    icablinkmetricsout = icablinkmetrics(EEG, 'ArtifactChannel', data_eog_blink, 'Alpha', 0.001, 'VisualizeData', 'False');
    if any(icablinkmetricsout.identifiedcomponents > 0)
        evidence.blinkmetrics(icablinkmetricsout.identifiedcomponents) = true;
    end
catch
    fprintf('The method has failed. Skipping...\n');
    icablinkmetricsout.identifiedcomponents = [];
    icablinkmetricsout.metrics = struct('corr_Pvalue', [], 'conv_Pvalue', [], 'perc_Pvalue', []);
end

EEG.ALSUTRECHT.ica.blink.icablinkmetrics.pval = single([icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]');
EEG.ALSUTRECHT.ica.blink.icablinkmetrics.bics = icablinkmetricsout.identifiedcomponents(:);
end

function [EEG, evidence] = gather_spatial_characteristics(EEG, cfg, evidence)
% GATHER_SPATIAL_CHARACTERISTICS Evaluates spatial smoothness/kurtosis across IC topographies.

if nargin < 2 || isempty(cfg), cfg = struct(); end
if nargin < 3 || isempty(evidence), evidence = struct(); end

fprintf('\n--------------------------------\n');
fprintf('Channel: Spatial characteristics\n');
fprintf('--------------------------------\n');

% threshold_spatial = 20; % kurtosis
threshold_spatial = 0.20; % focal

[spatial_smoothness, bad_ic] = estimate_spatialsmoothnes(EEG, threshold_spatial);

num_ica = size(EEG.icaact, 1);
evidence.spatial = false(num_ica, 1);
evidence.spatial(bad_ic) = true;

EEG.ALSUTRECHT.ica.channel.spatial_smoothness.focal = spatial_smoothness;
EEG.ALSUTRECHT.ica.channel.spatial_smoothness.bics  = bad_ic(:);

% -------------------------------------------------------------------------
% Diagnostic Figure: Metric Bar Plot with Threshold
% -------------------------------------------------------------------------
do_plot = ~isfield(cfg, 'do_plot') || cfg.do_plot;
if do_plot
    fh = figure('Color', 'w', 'Position', [150 150 950 360], 'Name', 'IC Spatial Characteristics');
    ax = axes(fh); hold(ax, 'on');

    ic_idx = 1:num_ica;
    vals   = spatial_smoothness(:);

    % Color clean components in muted gray, flagged components in red
    bar_colors = repmat([0.72 0.75 0.78], num_ica, 1);
    bar_colors(evidence.spatial, :) = repmat([0.85 0.20 0.15], sum(evidence.spatial), 1);

    b = bar(ax, ic_idx, vals, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.7);
    b.CData = bar_colors;

    % Threshold marker
    yline(ax, threshold_spatial, 'r--', ...
        sprintf('Threshold (%.0f)', threshold_spatial), ...
        'LineWidth', 1.3, 'FontSize', 10, 'LabelHorizontalAlignment', 'left');

    xlabel(ax, 'Independent Component (IC)', 'FontWeight', 'bold');
    ylabel(ax, 'Spatial Kurtosis / Smoothness', 'FontWeight', 'bold');
    title(ax, sprintf('IC Spatial Smoothness (%d / %d Flagged)', sum(evidence.spatial), num_ica), ...
        'FontSize', 11, 'FontWeight', 'bold');

    xlim(ax, [0.5, num_ica + 0.5]);
    box(ax, 'off');
    grid(ax, 'on');
    xticks(ax, 1:num_ica); xtickangle(ax, 90);

    % Auto-save if subject figures path is configured
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_channel_spatial'], [22 9]);
end

end

function [EEG, evidence] = gather_template_matches(EEG, templates_ica, evidence)
fprintf('\n--------------------------------\n');
fprintf('Template Matching (Blink, Saccade, Heart)\n');
fprintf('--------------------------------\n');

num_ica = size(EEG.icaact, 1);
evidence.template_blink   = false(num_ica, 1);
evidence.template_saccade = false(num_ica, 1);
evidence.template_heart   = false(num_ica, 1);

% Blink
cfg_base = struct('check_distribution', false, 'sim_thresh', 0.85, 'biosemi_blink', 'blink');
[bad_ic_blink, ~, ~, rep_blink] = match_ica_template(EEG.icawinv, EEG.chanlocs, templates_ica, 'blink', cfg_base);
evidence.template_blink(bad_ic_blink) = true;
EEG.ALSUTRECHT.ica.blink.TemplateCorr.bics   = bad_ic_blink(:);
EEG.ALSUTRECHT.ica.blink.TemplateCorr.report = rep_blink;

% Saccade
cfg_base = struct('check_distribution', false, 'sim_thresh', 0.85, 'biosemi_saccade', 'saccade');
[bad_ic_saccade, ~, ~, rep_saccade] = match_ica_template(EEG.icawinv, EEG.chanlocs, templates_ica, 'saccade', cfg_base);
evidence.template_saccade(bad_ic_saccade) = true;
EEG.ALSUTRECHT.ica.saccade.TemplateCorr.bics   = bad_ic_saccade(:);
EEG.ALSUTRECHT.ica.saccade.TemplateCorr.report = rep_saccade;

% Heart
cfg_base = struct('check_distribution', false, 'sim_thresh', 0.85);
[bad_ic_heart, ~, ~, rep_heart] = match_ica_template(EEG.icawinv, EEG.chanlocs, templates_ica, 'heart', cfg_base);
evidence.template_heart(bad_ic_heart) = true;
EEG.ALSUTRECHT.ica.heart.TemplateCorr.bics   = bad_ic_heart(:);
EEG.ALSUTRECHT.ica.heart.TemplateCorr.report = rep_heart;
end

function [EEG, evidence] = gather_spectral_metrics(EEG, cfg, evidence)
% GATHER_SPECTRAL_METRICS Evaluates atypical frequency roll-off / negative spectral
% integral across all ICs, flags bad components exceeding threshold, and renders
% a two-panel diagnostic display (Power Spectra + Decision Metric Space) along
% with bottom-aligned topoplots.

fprintf('\n--------------------------------\n');
fprintf('Bad (general): Spectral power\n');
fprintf('--------------------------------\n');

interest_range = [1 70];
censor_range   = [2 35];
threshold_power_prop    = 0.35;
threshold_power_slope_1 = -0.7; % EMG
threshold_power_slope_2 = -1.5; % EOG

data_ica              = EEG.icaact;
num_ica               = size(data_ica, 1);
evidence.spectral_bad = false(num_ica, 1);

% Estimate power spectral density across all ICs
[pow, frefAll] = pwelch(data_ica', size(data_ica, 2), [], size(data_ica, 2), EEG.srate);
freq = 1:0.5:100;
fftBins = zeros(size(pow, 2), size(freq, 2));

for i_freq = 1:length(freq)
    [~, index1] = min(abs(frefAll - (freq(i_freq) - 0.25)));
    [~, index2] = min(abs(frefAll - (freq(i_freq) + 0.25)));
    fftBins(:, i_freq) = mean(pow(index1:index2, :), 1)';
end

% Spectral regression
% [neg_integral, neg_prop, slope] = estimate_spectral_metrics(fftBins, freq, [1 70], [2 40]);
[neg_integral, neg_prop, slope] = estimate_spectral_metrics(fftBins, freq, interest_range, censor_range);

bad_ic_1 = neg_prop > threshold_power_prop;
bad_ic_2 = slope    > threshold_power_slope_1;
bad_ic_3 = slope    < threshold_power_slope_2;

bad_ic   = (bad_ic_1(:) & bad_ic_2(:)) | bad_ic_3(:);
bad_ic   = find(bad_ic);
evidence.spectral_bad(bad_ic) = true;

EEG.ALSUTRECHT.ica.bad.powerspectra.neg_integral = neg_integral;
EEG.ALSUTRECHT.ica.bad.powerspectra.neg_prop     = neg_prop;
EEG.ALSUTRECHT.ica.bad.powerspectra.slope        = slope;
EEG.ALSUTRECHT.ica.bad.powerspectra.bics         = bad_ic(:);

% -------------------------------------------------------------------------
% Visualisation Setup
% -------------------------------------------------------------------------
if ~cfg.figure.plot
    return;
end

n_flagged = length(bad_ic);
winv      = EEG.icawinv;
psd_db    = 10 * log10(fftBins);

if n_flagged <= 8
    colours = brewermap(max(3, n_flagged), 'Dark2');
elseif n_flagged <= 12
    colours = brewermap(n_flagged, 'Paired');
else
    colours = lines(n_flagged);
end

% -------------------------------------------------------------------------
% Dynamic Layout Configuration
% -------------------------------------------------------------------------
if n_flagged == 0
    % Compact single-row figure with 2 diagnostic panels side-by-side
    fig_width  = 1000;
    fig_height = 440;
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.figure.visible);

    t_top = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    ax_spec = nexttile(t_top);
    ax_meas = nexttile(t_top);

    plotX = 26;
    plotY = 12;
else
    max_cols_single_row = 8;
    if n_flagged <= max_cols_single_row
        n_cols      = n_flagged;
        n_topo_rows = 1;
    else
        n_cols      = 6;
        n_topo_rows = ceil(n_flagged / n_cols);
    end

    fig_width  = max(1000, 210 * n_cols + 100);
    fig_height = 500 + (210 * n_topo_rows);
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.figure.visible);

    % Master 2-row layout: Row 1 = Diagnostics; Row 2 = Topoplots
    t_main = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

    % Row 1 Sub-layout: 2 side-by-side diagnostic plots
    t_top = tiledlayout(t_main, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_top.Layout.Tile = 1;
    ax_spec = nexttile(t_top);
    ax_meas = nexttile(t_top);

    plotX = max(26, 3.8 * n_cols + 4);
    plotY = 10 + (4.5 * n_topo_rows);
end

% -------------------------------------------------------------------------
% Plot 1: Power Spectra (1-70 Hz)
% -------------------------------------------------------------------------
hold(ax_spec, 'on');

% Background: All IC spectra
h_all = plot(ax_spec, freq, psd_db', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);

% Flagged components
h_flagged   = gobjects(n_flagged, 1);
leg_entries = {'All ICs'};

for i = 1:n_flagged
    ic_idx = bad_ic(i);
    h_flagged(i) = plot(ax_spec, freq, psd_db(ic_idx, :), ...
        'Color', colours(i, :), 'LineWidth', 2.0);
    leg_entries{end+1} = sprintf('IC%d', ic_idx);
end

grid(ax_spec, 'on');
set(ax_spec, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
xlim(ax_spec, [1, 70]);
xlabel(ax_spec, 'Frequency (Hz)', 'FontWeight', 'bold');
ylabel(ax_spec, 'Power (10\cdotlog_{10} \muV^2/Hz)', 'FontWeight', 'bold');
title(ax_spec, 'Component Power Spectra', 'FontSize', 12, 'FontWeight', 'bold');

if n_flagged == 0
    subtitle(ax_spec, 'No ICs exceeded spectral criteria', 'FontAngle', 'italic', 'Color', [0.45, 0.45, 0.45]);
    legend(ax_spec, h_all(1), 'All ICs', 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
else
    n_leg_cols = min(3, ceil(length(leg_entries) / 8));
    legend(ax_spec, [h_all(1); h_flagged], leg_entries, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 8, 'NumColumns', n_leg_cols);
end
hold(ax_spec, 'off');

% -------------------------------------------------------------------------
% Plot 2: Decision Metric Space (Proportion vs Slope)
% -------------------------------------------------------------------------
hold(ax_meas, 'on');

% Axis ranges with padding
% x_min = max(0, min(neg_prop) - 0.05);
x_min = 0.05;
x_max = max(threshold_power_prop + 0.15, max(neg_prop) + 0.05);
y_min = min(threshold_power_slope_1 - 0.4, min(slope) - 0.1);
y_max = max(threshold_power_slope_1 + 0.3, max(slope) + 0.1);

% Shaded Rejection Quadrant (Upper-Right)
patch(ax_meas, [threshold_power_prop, x_max, x_max, threshold_power_prop], ...
    [threshold_power_slope_1, threshold_power_slope_1, y_max, y_max], ...
    [1.0, 0.92, 0.92], 'EdgeColor', 'none', 'DisplayName', 'Rejection Zone');

% Decision threshold lines
xline(ax_meas, threshold_power_prop, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2, 'DisplayName', 'Threshold 1');
yline(ax_meas, threshold_power_slope_1, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2, 'DisplayName', 'Threshold 2');

% Clean ICs (unflagged)
clean_mask = ~evidence.spectral_bad;
clean_idx  = find(clean_mask);

scatter(ax_meas, neg_prop(clean_mask), slope(clean_mask), 35, ...
    [0.65, 0.65, 0.65], 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Clean ICs');

% Vectorised grey labels for clean ICs to audit borderline / false negatives
clean_labels = cellstr(string(clean_idx));
text(ax_meas, neg_prop(clean_mask) + 0.005, slope(clean_mask), clean_labels, ...
    'Color', [0.55, 0.55, 0.55], 'FontSize', 7.5, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'HandleVisibility', 'off');

% Flagged ICs (colour-matched with bold labels)
for i = 1:n_flagged
    ic_idx = bad_ic(i);
    scatter(ax_meas, neg_prop(ic_idx), slope(ic_idx), 75, colours(i, :), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    text(ax_meas, neg_prop(ic_idx) + 0.008, slope(ic_idx), sprintf('%d', ic_idx), ...
        'Color', colours(i, :), 'FontWeight', 'bold', 'FontSize', 9.5, ...
        'VerticalAlignment', 'middle');
end

xlim(ax_meas, [x_min, x_max]);
ylim(ax_meas, [y_min, y_max]);
grid(ax_meas, 'on');
set(ax_meas, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
xlabel(ax_meas, 'Negative Spectral Proportion (prop)', 'FontWeight', 'bold');
ylabel(ax_meas, 'Power Slope (slope)', 'FontWeight', 'bold');
title(ax_meas, 'Spectral Decision Space', 'FontSize', 12, 'FontWeight', 'bold');
subtitle(ax_meas, sprintf('Flagged if prop > %.2f & slope > %.2f', threshold_power_prop, threshold_power_slope_1), ...
    'FontAngle', 'italic', 'Color', [0.45, 0.45, 0.45]);

legend(ax_meas, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
hold(ax_meas, 'off');
% -------------------------------------------------------------------------
% Panel 3: Topoplots of Flagged Components (Full Bottom Width)
% -------------------------------------------------------------------------
if n_flagged > 0
    t_bot = tiledlayout(t_main, n_topo_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_bot.Layout.Tile = 2;

    for i = 1:n_flagged
        ax_topo = nexttile(t_bot);
        ic_idx     = bad_ic(i);
        pc_weights = winv(:, ic_idx);
        cmax       = max(abs(pc_weights));
        if cmax == 0, cmax = 1; end

        title_str = sprintf('IC%d\n(%.2f, %.2f)', ic_idx, neg_prop(ic_idx), slope(ic_idx));
        mytopoplot(pc_weights, [], '', ax_topo, [-cmax, cmax]);
        title(ax_topo, title_str, 'Color', colours(i, :), 'FontWeight', 'bold', 'FontSize', 8);
    end
end

% -------------------------------------------------------------------------
% Save Figure
% -------------------------------------------------------------------------
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_spectral_bad'], [plotX, plotY]);
end

function [EEG, evidence] = gather_temporal_metrics(EEG, cfg, evidence)
% GATHER_TEMPORAL_METRICS Evaluates temporal activation properties across all
% ICs (sample-to-sample high-frequency roughness and kurtosis), flags artifactual
% components exceeding thresholds, and renders a diagnostic metric space plot
% (with kurtosis display capped to preserve scale) along with bottom-aligned topoplots.

fprintf('\n--------------------------------\n');
fprintf('Bad (general): Temporal activation properties\n');
fprintf('--------------------------------\n');

% -------------------------------------------------------------------------
% 1. Parameter Validation & Configuration Defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg, 'figure'),               cfg.figure = struct(); end
if ~isfield(cfg.figure, 'plot'),          cfg.figure.plot = true; end
if ~isfield(cfg.figure, 'visible'),       cfg.figure.visible = 'on'; end

if ~isfield(cfg, 'kurtosis_thresh'),      cfg.kurtosis_thresh = 5.25;  end
if ~isfield(cfg, 'hf_rough_thresh'),      cfg.hf_rough_thresh = 0.35;  end
if ~isfield(cfg, 'kurtosis_cap'),         cfg.kurtosis_cap    = 40;    end % Display cap
if ~isfield(cfg, 'temporal_mode'),        cfg.temporal_mode   = 'or';  end % 'or' | 'and'

% Ensure ICA activations are computed
if isempty(EEG.icaact)
    ch_idx = EEG.icachansind;
    EEG.icaact = (EEG.icaweights * EEG.icasphere) * reshape(EEG.data(ch_idx, :, :), length(ch_idx), []);
end

data_ica = EEG.icaact;
num_ica  = size(data_ica, 1);
evidence.temporal_bad = false(num_ica, 1);

% Reshape epoched data to 2D continuous activations if necessary
if ndims(data_ica) == 3
    data_ica_2d = reshape(data_ica, num_ica, []);
else
    data_ica_2d = data_ica;
end

% -------------------------------------------------------------------------
% 2. Compute Temporal Metrics (Kurtosis & HF Roughness)
% -------------------------------------------------------------------------
kurt         = zeros(num_ica, 1);
hf_roughness = zeros(num_ica, 1);

for i_ic = 1:num_ica
    act_k = double(data_ica_2d(i_ic, :));

    % True kurtosis calculation
    kurt(i_ic) = kurtosis(act_k);

    % High-frequency roughness: Var(dx) / (2 * Var(x))
    diff_act = diff(act_k);
    var_act  = var(act_k);
    hf_roughness(i_ic) = var(diff_act) / (2 * var_act + eps);
end

% Decision gates evaluated on true uncapped values
bad_ic_1 = kurt > cfg.kurtosis_thresh;
bad_ic_2 = hf_roughness > cfg.hf_rough_thresh;

if strcmpi(cfg.temporal_mode, 'and')
    bad_mask = bad_ic_1(:) & bad_ic_2(:);
else
    bad_mask = bad_ic_1(:) | bad_ic_2(:);
end

bad_ic = find(bad_mask);
evidence.temporal_bad(bad_ic) = true;

% Log uncapped values to EEG struct
EEG.ALSUTRECHT.ica.bad.temporal.kurtosis     = kurt;
EEG.ALSUTRECHT.ica.bad.temporal.hf_roughness = hf_roughness;
EEG.ALSUTRECHT.ica.bad.temporal.bics         = bad_ic(:);

% -------------------------------------------------------------------------
% 3. Visualisation Setup
% -------------------------------------------------------------------------
if ~cfg.figure.plot
    return;
end

n_flagged = length(bad_ic);
winv      = EEG.icawinv;

if n_flagged <= 8
    colours = brewermap(max(3, n_flagged), 'Dark2');
elseif n_flagged <= 12
    colours = brewermap(n_flagged, 'Paired');
else
    colours = lines(n_flagged);
end

% Cap kurtosis strictly for plotting to prevent axis squashing
kurt_plot = min(kurt, cfg.kurtosis_cap);

% -------------------------------------------------------------------------
% Dynamic Layout Configuration
% -------------------------------------------------------------------------
if n_flagged == 0
    fig_width  = 580;
    fig_height = 450;
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.figure.visible);

    t_main  = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_meas = nexttile(t_main);

    plotX = 15;
    plotY = 11;
else
    max_cols_single_row = 8;
    if n_flagged <= max_cols_single_row
        n_cols      = n_flagged;
        n_topo_rows = 1;
    else
        n_cols      = 6;
        n_topo_rows = ceil(n_flagged / n_cols);
    end

    fig_width  = max(650, 210 * n_cols + 80);
    fig_height = 360 + (210 * n_topo_rows);
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.figure.visible);

    t_main = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    ax_meas = nexttile(t_main, 1);

    plotX = max(16, 3.8 * n_cols + 3);
    plotY = 9.5 + (4.5 * n_topo_rows);
end

% -------------------------------------------------------------------------
% Panel 1: Decision Metric Space (Roughness vs Kurtosis)
% -------------------------------------------------------------------------
hold(ax_meas, 'on');

% Axis ranges with padding
x_min = max(0, min(hf_roughness) - 0.05);
x_max = max(cfg.hf_rough_thresh + 0.20, max(hf_roughness) + 0.05);
y_min = max(0, min(kurt_plot) - 0.8);
y_max = cfg.kurtosis_cap + 5.0;

% Shaded Rejection Zone
if strcmpi(cfg.temporal_mode, 'and')
    patch(ax_meas, [cfg.hf_rough_thresh, x_max, x_max, cfg.hf_rough_thresh], ...
        [cfg.kurtosis_thresh, cfg.kurtosis_thresh, y_max, y_max], ...
        [1.0, 0.92, 0.92], 'EdgeColor', 'none', 'DisplayName', 'Rejection Zone');
else
    % L-shaped rejection region for OR logic
    patch(ax_meas, [cfg.hf_rough_thresh, x_max, x_max, cfg.hf_rough_thresh], ...
        [y_min, y_min, y_max, y_max], ...
        [1.0, 0.92, 0.92], 'EdgeColor', 'none', 'DisplayName', 'Rejection Zone');
    patch(ax_meas, [x_min, cfg.hf_rough_thresh, cfg.hf_rough_thresh, x_min], ...
        [cfg.kurtosis_thresh, cfg.kurtosis_thresh, y_max, y_max], ...
        [1.0, 0.92, 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Decision threshold lines
xline(ax_meas, cfg.hf_rough_thresh, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2, 'DisplayName', 'Roughness Threshold');
yline(ax_meas, cfg.kurtosis_thresh, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2, 'DisplayName', 'Kurtosis Threshold');

% Clean ICs (unflagged)
clean_mask = ~evidence.temporal_bad;
clean_idx  = find(clean_mask);

scatter(ax_meas, hf_roughness(clean_mask), kurt_plot(clean_mask), 35, ...
    [0.65, 0.65, 0.65], 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Clean ICs');

% Vectorised grey labels for clean ICs to audit false negatives
clean_labels = cellstr(string(clean_idx));
text(ax_meas, hf_roughness(clean_mask) + 0.005, kurt_plot(clean_mask), clean_labels, ...
    'Color', [0.55, 0.55, 0.55], 'FontSize', 7.5, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'HandleVisibility', 'off');

% Flagged ICs (colour-matched with bold labels)
for i = 1:n_flagged
    ic_idx = bad_ic(i);
    scatter(ax_meas, hf_roughness(ic_idx), kurt_plot(ic_idx), 75, colours(i, :), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Append an arrow or indicator if kurtosis was clipped by the display cap
    if kurt(ic_idx) > cfg.kurtosis_cap
        lbl_str = sprintf('%d (\\uparrow)', ic_idx);
    else
        lbl_str = sprintf('%d', ic_idx);
    end

    text(ax_meas, hf_roughness(ic_idx) + 0.008, kurt_plot(ic_idx), lbl_str, ...
        'Color', colours(i, :), 'FontWeight', 'bold', 'FontSize', 9.5, ...
        'VerticalAlignment', 'middle');
end

xlim(ax_meas, [x_min, x_max]);
ylim(ax_meas, [y_min, y_max]);
grid(ax_meas, 'on');
set(ax_meas, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
xlabel(ax_meas, 'HF Roughness (Var(\Delta x) / 2\cdotVar(x))', 'FontWeight', 'bold');
ylabel(ax_meas, sprintf('Kurtosis (Capped at %g)', cfg.kurtosis_cap), 'FontWeight', 'bold');
title(ax_meas, 'Temporal Decision Space', 'FontSize', 12, 'FontWeight', 'bold');
subtitle(ax_meas, sprintf('Flagged if Kurtosis > %.1f %s Roughness > %.2f', ...
    cfg.kurtosis_thresh, upper(cfg.temporal_mode), cfg.hf_rough_thresh), ...
    'FontAngle', 'italic', 'Color', [0.45, 0.45, 0.45]);

legend(ax_meas, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
hold(ax_meas, 'off');

% -------------------------------------------------------------------------
% Panel 2: Topoplots of Flagged Components (Full Bottom Width)
% -------------------------------------------------------------------------
if n_flagged > 0
    t_bot = tiledlayout(t_main, n_topo_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_bot.Layout.Tile = 2;

    for i = 1:n_flagged
        ax_topo = nexttile(t_bot);
        ic_idx     = bad_ic(i);
        pc_weights = winv(:, ic_idx);
        cmax       = max(abs(pc_weights));
        if cmax == 0, cmax = 1; end

        % Displays the true uncapped kurtosis on the individual topoplots
        title_str = sprintf('IC%d\n(K=%.1f, R=%.2f)', ic_idx, kurt(ic_idx), hf_roughness(ic_idx));
        mytopoplot(pc_weights, [], '', ax_topo, [-cmax, cmax]);
        title(ax_topo, title_str, 'Color', colours(i, :), 'FontWeight', 'bold', 'FontSize', 8);
    end
end

% -------------------------------------------------------------------------
% Save Figure
% -------------------------------------------------------------------------
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_temporal_bad'], [plotX, plotY]);

end

%% ========================================================================
% SYNTHESIS AND REPORTING
% ========================================================================

function EEG = synthesise_evidence(EEG, evidence, templates_ica, cfg)
fprintf('\n--------------------------------\n');
fprintf('Combining and Synthesising Detected ICs\n');
fprintf('--------------------------------\n');

data_ica = EEG.icaact;
num_ica  = size(data_ica, 1);
icawinv  = EEG.icawinv;
ICLabel_struct = EEG.ALSUTRECHT.ica.ICLabel;

% 1. Blinks
blink_votes        = evidence.corr_veog + evidence.blinkmetrics + evidence.template_blink;
candidate_blinks   = find(blink_votes >= 1);
candidate_blinks   = is_false_veog(candidate_blinks, icawinv, templates_ica, ICLabel_struct);
ICsMostLikelyBlink = false(num_ica, 1);
ICsMostLikelyBlink(candidate_blinks) = true;

% 2. Saccades
saccade_votes        = evidence.corr_heog | evidence.template_saccade;
candidate_saccades   = find(saccade_votes);
candidate_saccades   = is_false_heog(candidate_saccades, icawinv, templates_ica, ICLabel_struct);
ICsMostLikelySaccade = false(num_ica, 1);
ICsMostLikelySaccade(candidate_saccades) = true;

% Eye aggregate
ICsMostLikelyEye = ICsMostLikelyBlink | ICsMostLikelySaccade | evidence.iclabel_eye;

% 3. Muscle
muscle_votes        = evidence.iclabel_muscle | evidence.emg_slope;
candidate_muscle    = find(muscle_votes);
candidate_muscle    = is_false_emg(candidate_muscle, icawinv, templates_ica, ICLabel_struct);
ICsMostLikelyMuscle = false(num_ica, 1);
ICsMostLikelyMuscle(candidate_muscle) = true;

% Complex Overlaps
ICsMostLikelyComplex = ICsMostLikelyMuscle & ICsMostLikelyEye;
ICsMostLikelyMuscle(ICsMostLikelyComplex) = false;
ICsMostLikelyEye(ICsMostLikelyComplex)    = false;

% 4. Channel Noise
channel_votes = evidence.iclabel_channel | evidence.spatial;
ICsMostLikelyChannel = channel_votes;
ICsMostLikelyChannel(ICsMostLikelyEye | ICsMostLikelyMuscle | ICsMostLikelyComplex) = false;

% 5. Heart
ICsMostLikelyHeart = detect_heart_synthesised(EEG, evidence, templates_ica, cfg);

% Enforce mutual exclusivity
ICsMostLikelyEye(ICsMostLikelyHeart)     = false;
ICsMostLikelyMuscle(ICsMostLikelyHeart)  = false;
ICsMostLikelyComplex(ICsMostLikelyHeart) = false;
ICsMostLikelyChannel(ICsMostLikelyHeart) = false;

% 6. General Bad
% candidate_bad            = find(evidence.spectral_bad | evidence.temporal_bad);
candidate_bad            = find(evidence.spectral_bad);
false_bad                = is_likely_brain(candidate_bad, ICLabel_struct, 0.60);
candidate_bad(false_bad) = [];
ICsMostLikelyBad = false(num_ica, 1);
ICsMostLikelyBad(candidate_bad) = true;

% Write final structures
EEG.ALSUTRECHT.ica.final.eye     = ICsMostLikelyEye;
EEG.ALSUTRECHT.ica.final.muscle  = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.final.complex = ICsMostLikelyComplex;
EEG.ALSUTRECHT.ica.final.channel = ICsMostLikelyChannel;
EEG.ALSUTRECHT.ica.final.heart   = ICsMostLikelyHeart;
EEG.ALSUTRECHT.ica.final.genbad  = ICsMostLikelyBad;

% Summary classification mapping
% 1: Brain, 2: Muscle, 3: Eye, 4: Heart, 6: Channel, 7: Other/Bad
EEG.ALSUTRECHT.ica.final.report                       = ICLabel_struct.cvec; % Basis
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyBad)     = 7;                   % General marker of artifact components
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyMuscle)  = 2;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyEye)     = 3;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyHeart)   = 4;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyChannel) = 6;

end

function EEG = compute_vaf_report(EEG, evidence, cfg)
% COMPUTE_VAF_REPORT Computes individual and group variance accounted for (VAF),
% logs results to EEG.ALSUTRECHT.ica.final.var, and renders the combined ICA
% evidence and classification dashboard figure.

num_ica = size(EEG.icaact, 1);
num_ica_relevant = min(25, num_ica);
fprintf('\nVariance Accounted For (VAF) report for first %d ICs...\n', num_ica_relevant);

% Extract continuous 2D channel data
if isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    ch_idx = EEG.icachansind;
elseif isfield(EEG.ALSUTRECHT.ica, 'icachansind')
    ch_idx = EEG.ALSUTRECHT.ica.icachansind;
else
    ch_idx = 1:size(EEG.data, 1);
end

ch_data = reshape(EEG.data(ch_idx, :, :), length(ch_idx), []);
icaact  = reshape(EEG.icaact, num_ica, []);

% Single IC VAF across the first K components
vaf_per_ic = zeros(1, num_ica_relevant);
for i_ic = 1:num_ica_relevant
    vaf_per_ic(i_ic) = get_vaf(ch_data, icaact, EEG.icawinv, i_ic);
end

% Cumulative VAF for relevant components
mask_relevant = false(num_ica, 1);
mask_relevant(1:num_ica_relevant) = true;
var_relevant = get_vaf(ch_data, icaact, EEG.icawinv, find(mask_relevant));
fprintf('Within the first %d ICs (scalp variance accounted for = %.2f%%):\n', num_ica_relevant, var_relevant);

% Class breakdown within relevant window
cat_labels = {'brain', 'muscle', 'eye', 'heart', 'line', 'channel', 'other'};
report_rel = EEG.ALSUTRECHT.ica.final.report(1:num_ica_relevant);
I = EEG.ALSUTRECHT.ica.final.report(:)';

for i_cat = 1:length(cat_labels)
    pct_comp = mean(report_rel == i_cat) * 100;
    vaf_val  = get_vaf(ch_data, icaact, EEG.icawinv, find(I == i_cat));
    EEG.ALSUTRECHT.ica.final.var.(cat_labels{i_cat}) = vaf_val;
    fprintf('%-10s components: %2.0f%% (scalp variance accounted for = %.2f%%)\n', cat_labels{i_cat}, round(pct_comp), vaf_val);
end

% General bad components
vaf_bad = get_vaf(ch_data, icaact, EEG.icawinv, find(EEG.ALSUTRECHT.ica.final.genbad(:)'));
EEG.ALSUTRECHT.ica.final.var.genbad = vaf_bad;
pct_bad_rel = mean(EEG.ALSUTRECHT.ica.final.genbad(1:num_ica_relevant)) * 100;
fprintf('%-10s components: %2.0f%% (scalp variance accounted for = %.2f%%)\n', 'bad', round(pct_bad_rel), vaf_bad);

% Render Unified ICA Evidence Dashboard
if cfg.figure.plot
    render_ica_summary_dashboard(EEG, evidence, num_ica_relevant, vaf_per_ic, cfg);
end

end

function render_ica_summary_dashboard(EEG, evidence, num_ica_relevant, vaf_per_ic, cfg)
% Renders an executive 3-panel layout capped at the first 60 ICs:
%   Panel 1: Evidence Matrix Heatmap across all gathering modules (with x-ticks)
%   Panel 2: Final Synthesised Class Strip with discrete colour legend
%   Panel 3: VAF Bar Chart for top components, colour-coded by verdict

% 1. Enforce strict 60-component window
num_ica  = size(EEG.icaweights, 1);
num_disp = min(60, num_ica);

% 2. Figure initialisation with wider aspect ratio
fh = figure('Color', 'w', 'Position', [50, 50, 1650, 750], ...
    'Visible', cfg.figure.visible, 'Name', 'ICA Verification Dashboard');
t = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Prepare Evidence Heatmap Matrix (16 metrics x n_disp ICs)
ev_fields = {'iclabel_eye', 'iclabel_muscle', 'iclabel_heart', 'iclabel_channel', ...
    'corr_veog', 'corr_heog', 'corr_ecg', 'ecg_ctps', 'ecg_erp', ...
    'emg_slope', 'blinkmetrics', 'spatial', 'template_blink', ...
    'template_saccade', 'template_heart', 'spectral_bad'};

ev_labels = {'ICL Eye', 'ICL Musc', 'ICL Heart', 'ICL Chan', ...
    'Corr VEOG', 'Corr HEOG', 'Corr ECG', 'CTPS pK', 'ECG ERP', ...
    'EMG Slope', 'BlinkMetr', 'Spatial Kurt', 'Tmpl Blink', ...
    'Tmpl Sacc', 'Tmpl Heart', 'Spectral Bad'};

ev_mat = zeros(length(ev_fields), num_disp);
for f = 1:length(ev_fields)
    if isfield(evidence, ev_fields{f})
        vec = evidence.(ev_fields{f})(:)';
        ev_mat(f, :) = double(vec(1:num_disp));
    end
end

% Category palette definitions [R G B]
% 1: Brain, 2: Muscle, 3: Eye, 4: Heart, 5: Line, 6: Channel, 7: Bad/Other
class_colors = [ ...
    0.20, 0.65, 0.35; ... % Brain (Green)
    0.85, 0.45, 0.15; ... % Muscle (Orange)
    0.20, 0.45, 0.75; ... % Eye (Blue)
    0.75, 0.15, 0.20; ... % Heart (Red)
    0.50, 0.50, 0.50; ... % Line (Grey)
    0.55, 0.25, 0.65; ... % Channel (Purple)
    0.30, 0.30, 0.30  ... % Other/Bad (Dark Grey)
    ];
legend_labels = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other / Bad'};

% -------------------------------------------------------------------------
% PANEL 1: Evidence Matrix (Top 2 tiles) with X-ticks
% -------------------------------------------------------------------------
ax1 = nexttile(t, [2 1]);
imagesc(ax1, 1:num_disp, 1:length(ev_fields), ev_mat);
colormap(ax1, [0.94, 0.94, 0.94; 0.15, 0.35, 0.55]); % Light Grey (Pass) vs Deep Slate (Flagged)
set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 8, 'YDir', 'normal');
yticks(ax1, 1:length(ev_fields));
yticklabels(ax1, ev_labels);

% Show full component indices along the top edge
set(ax1, 'XAxisLocation', 'top');
xticks(ax1, 1:num_disp);
xticklabels(ax1, 1:num_disp);
set(ax1, 'TickDir', 'out');

title(ax1, sprintf('%s: Independent Evidence Tracker (ICs 1–%d)', EEG.ALSUTRECHT.subject.id, num_disp), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid(ax1, 'on');
set(ax1, 'GridAlpha', 0.2, 'GridLineStyle', '-');
xlim(ax1, [0.5, num_disp + 0.5]);
axis tight;

% -------------------------------------------------------------------------
% PANEL 2: Synthesised Class Strip (1 tile) with Discrete Swatch Legend
% -------------------------------------------------------------------------
ax2 = nexttile(t, [1 1]);
final_classes = EEG.ALSUTRECHT.ica.final.report(1:num_disp)';
imagesc(ax2, 1:num_disp, 1, final_classes);
colormap(ax2, class_colors);
clim(ax2, [1, 7]);
set(ax2, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 8, 'YDir', 'normal');
yticks(ax2, 1);
yticklabels(ax2, {'Verdict'});
xticks(ax2, 1:num_disp);
xticklabels(ax2, 1:num_disp);
set(ax2, 'TickDir', 'out');
% xlabel(ax2, sprintf('Independent Components (ICs 1–%d)', n_disp), 'FontSize', 10, 'FontWeight', 'bold');
% title(ax2, 'Synthesised Classification Verdict', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax2, [0.5, num_disp + 0.5]);

% Create dummy swatch markers for horizontal legend
hold(ax2, 'on');
h_leg = gobjects(length(legend_labels), 1);
for i_c = 1:length(legend_labels)
    h_leg(i_c) = plot(ax2, NaN, NaN, 's', ...
        'MarkerFaceColor', class_colors(i_c, :), ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25], ...
        'MarkerSize', 8, 'LineWidth', 0.5);
end
hold(ax2, 'off');

legend(ax2, h_leg, legend_labels, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'Box', 'on', 'FontSize', 9);

% -------------------------------------------------------------------------
% PANEL 3: Top VAF Breakdown (1 tile)
% -------------------------------------------------------------------------
ax3 = nexttile(t, [1 1]);
n_vaf_bars = min(num_ica_relevant, num_disp);
b = bar(ax3, 1:n_vaf_bars, vaf_per_ic(1:n_vaf_bars), 'BarWidth', 0.6, 'EdgeColor', 'none');
b.FaceColor = 'flat';
for i_b = 1:n_vaf_bars
    cls_idx = max(1, min(7, final_classes(i_b)));
    b.CData(i_b, :) = class_colors(cls_idx, :);
end
grid(ax3, 'on');
set(ax3, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 8, 'GridAlpha', 0.4);
ylabel(ax3, 'VAF (%)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(ax3, 'Component', 'FontSize', 10, 'FontWeight', 'bold');
title(ax3, sprintf('Individual Scalp Variance Accounted For (Cumulative Top %d = %.1f%%)', ...
    n_vaf_bars, sum(vaf_per_ic(1:n_vaf_bars))), 'FontSize', 10, 'FontWeight', 'bold');
xlim(ax3, [0.5, n_vaf_bars + 0.5]);
xticks(ax3, 1:n_vaf_bars);

linkaxes([ax1, ax2], 'x');

% 3. Save wide figure (36 cm x 16 cm)
plotX = 36;
plotY = 16;
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_summary_dashboard'], [plotX, plotY]);
end

%% ========================================================================
% VALIDATION HELPERS & LOCAL ROUTINES
% ========================================================================
function vaf = get_vaf(ch_data, icaact, icawinv, comp_idx)
comp_idx = comp_idx(:)';
if isempty(comp_idx)
    vaf = 0;
    return;
end

if exist('compvar', 'file') == 2
    [~, vaf] = compvar(ch_data, icaact, icawinv, comp_idx);
else
    tot_var  = sum(var(ch_data, 0, 2));
    proj     = icawinv(:, comp_idx) * icaact(comp_idx, :);
    resid_var = sum(var(ch_data - proj, 0, 2));
    vaf      = max(0, (1 - (resid_var / tot_var)) * 100);
end
end

function bad_ic = is_false_veog(bad_ic, ICs, templates_ica, ICLabel)
if isempty(bad_ic); return; end

blink_templates = [templates_ica.Blinkweights0, templates_ica.Blinkweights1];
max_sim_blink = calc_max_sim(ICs(:, bad_ic), blink_templates, 'cosine');

false_no_template = max_sim_blink < 0.80;
false_is_brain    = is_likely_brain(bad_ic, ICLabel, 0.60);
false_all         = false_no_template(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[VEOG Validation] Pruned %d/%d false-positive VEOG ICs.\n', sum(false_all), length(false_all));
end

function ICsMostLikelyHeart = detect_heart_synthesised(EEG, evidence, templates_ica, cfg)
% DETECT_HEART_SYNTHESISED Synthesises ICLabel, spatial template matching,
% and temporal ECG metrics (R-peak ERP, CTPS, cross-correlation) to identify
% both primary far-field cardiac dipoles and peripheral rim/neck split components.

is_ecg_available = strcmpi(EEG.ALSUTRECHT.subject.ecg, 'recorded');

data_ica = EEG.icaact;
num_ica  = size(data_ica, 1);
ICsMostLikelyHeart = false(num_ica, 1);

iclabel_class = EEG.ALSUTRECHT.ica.ICLabel.cvec;
iclabel_probs = EEG.ALSUTRECHT.ica.ICLabel.pvec;

% Check spectral match
cfg_match = struct('sim_thresh', 0.55, 'metric', 'cosine');
[~, sim_scores] = match_ica_template(EEG.icawinv, EEG.chanlocs, templates_ica, 'heart', cfg_match);

% Check how many candidates there are
cand_ecg_ics = ICsMostLikelyHeart | evidence.ecg_erp | evidence.ecg_ctps | evidence.template_heart;

if any(cand_ecg_ics)
    % Check their spectra
    [has_alpha_peak_tmp, has_ecg_spectrum_tmp] = check_ecg_candidate_spectra(EEG, cand_ecg_ics, EEG.srate, cfg);

    has_alpha_peak   = false(size(cand_ecg_ics));
    has_ecg_spectrum = false(size(cand_ecg_ics));
    cand_ecg_ics_tmp = find(cand_ecg_ics);
    has_alpha_peak(cand_ecg_ics_tmp(has_alpha_peak_tmp))     = true;
    has_ecg_spectrum(cand_ecg_ics_tmp(has_ecg_spectrum_tmp)) = true;

    for i_comp = 1:num_ica
        % 1. Direct ICLabel Classification
        if iclabel_class(i_comp) == 4
            ICsMostLikelyHeart(i_comp) = true;
            continue;
        end

        if is_ecg_available
            % 2. Ground-Truth QRS Averaging (R-peak ERP)
            % R-peak triggered averaging across hundreds of epochs is definitive temporal
            % proof. Peripheral rim/neck components (split components like IC 15 and IC 18)
            % pick up local volume conduction and will fail global far-field template matching
            % (sim_scores < 0.55), but remain unequivocally cardiac.
            if cand_ecg_ics(i_comp) && ~has_alpha_peak(i_comp) && has_ecg_spectrum(i_comp)
                ICsMostLikelyHeart(i_comp) = true;
                continue;
            end

            % 3. Cycle-Triggered Probability / Continuous Correlation
            % These metrics carry higher risk of spurious correlation from slow drift,
            % so they still require spatial template confirmation.
            if evidence.ecg_ctps(i_comp) && sim_scores(i_comp) > 0.90
                ICsMostLikelyHeart(i_comp) = true;
            elseif evidence.corr_ecg(i_comp) && sim_scores(i_comp) > 0.90
                ICsMostLikelyHeart(i_comp) = true;
            end

        else
            % 4. Fallback: No ECG Lead Available (Spatial Topography & Far-Field Physics Only)
            is_strict_spatial = sim_scores(i_comp) > 0.90;

            if is_strict_spatial && ~has_alpha_peak(i_comp) && has_ecg_spectrum(i_comp)
                ICsMostLikelyHeart(i_comp) = true;
            end
        end


        % if is_ecg_available
        %     % 2. Ground-Truth QRS Averaging (R-peak ERP)
        %     % R-peak triggered averaging across hundreds of epochs is definitive temporal
        %     % proof. Peripheral rim/neck components (split components like IC 15 and IC 18)
        %     % pick up local volume conduction and will fail global far-field template matching
        %     % (sim_scores < 0.55), but remain unequivocally cardiac.
        %     if evidence.ecg_erp(ic) && evidence.ecg_ctps(ic)
        %         ICsMostLikelyHeart(ic) = true;
        %         continue;
        %     end
        %
        %     % 3. Cycle-Triggered Probability / Continuous Correlation
        %     % These metrics carry higher risk of spurious correlation from slow drift,
        %     % so they still require spatial template confirmation.
        %     if evidence.ecg_ctps(ic) && sim_scores(ic) > 0.80
        %         ICsMostLikelyHeart(ic) = true;
        %     elseif evidence.corr_ecg(ic) && sim_scores(ic) > 0.80
        %         ICsMostLikelyHeart(ic) = true;
        %     end
        %
        % else
        %     % 4. Fallback: No ECG Lead Available (Spatial Topography & Far-Field Physics Only)
        %     is_strict_spatial = (sim_scores(ic) >= 0.95);
        %     ic_weights        = EEG.icawinv(:, ic);
        %     focality          = max(abs(ic_weights)) / norm(ic_weights);
        %     is_far_field      = focality < 0.28;
        %
        %     if is_strict_spatial && is_far_field
        %         ICsMostLikelyHeart(ic) = true;
        %     end
        % end
    end

    % Note: Do not impose an arbitrary cap (e.g. max 2 components). High-density montages
    % regularly split the cardiac vector into 1 global dipole and 1 to 2 perimeter rim components.

end

end

function [has_alpha_peak, has_ecg_spectrum] = check_ecg_candidate_spectra(EEG, cand_ecg_ics, srate, cfg)
% CHECK_ECG_CANDIDATE_SPECTRA Plots power spectra of candidate ECG independent
% components against the background of all ICs to differentiate true cardiac
% activity (steep 1/f roll-off, power concentrated < 35 Hz) from myogenic
% chimeras (elevated broadband plateau between 40 and 100 Hz).
%
% Inputs:
%   EEG          - EEGLAB structure containing .icaact, .icawinv, and metadata
%   cand_ecg_ics - Vector of IC indices (numeric) or logical mask [num_ica x 1]
%   srate        - Sampling frequency in Hz
%   cfg          - (Optional) configuration struct:
%                    .max_freq     : Upper frequency display limit (default: 80 Hz)
%                    .win_sec      : Welch window length in seconds (default: 2 s)
%                    .plot_visible : 'on' (default) or 'off'
%                    .figure.plot  : true (default) or false

% -------------------------------------------------------------------------
% 1. Input Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 4 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg, 'max_freq'),     cfg.max_freq     = 60;    end
if ~isfield(cfg, 'win_sec'),      cfg.win_sec      = 8.0;   end
if ~isfield(cfg, 'plot_visible'), cfg.plot_visible = 'on';  end

if ~isfield(cfg, 'figure'), cfg.figure = struct(); end
if ~isfield(cfg.figure, 'plot')
    cfg.figure.plot = true;
end
if ~isfield(cfg.figure, 'visible')
    cfg.figure.visible = cfg.plot_visible;
end

if islogical(cand_ecg_ics)
    cand_ecg_ics = find(cand_ecg_ics);
end
cand_ecg_ics = cand_ecg_ics(:)';
n_cands  = length(cand_ecg_ics);

if n_cands == 0
    warning('No candidate ICs provided to plot.');
    has_alpha_peak   = NaN;
    has_ecg_spectrum = NaN;
    return;
end

% -------------------------------------------------------------------------
% 2. Estimate Power Spectral Density (Smoothed Welch)
% -------------------------------------------------------------------------
data_ica    = EEG.icaact;
win_samples = min(size(data_ica, 2), round(cfg.win_sec * srate));
n_overlap   = round(win_samples / 2);
n_fft       = win_samples;

[pow, freq] = pwelch(data_ica', win_samples, n_overlap, n_fft, srate);
psd_db      = 10 * log10(pow'); % [num_ica x freq_bins]

% Restrict to frequency range
f_mask = freq >= 1 & freq <= cfg.max_freq;
freq   = freq(f_mask);
psd_db = psd_db(:, f_mask);

% -------------------------------------------------------------------------
% 3. Check for Spectral Properties
% -------------------------------------------------------------------------
% Check for 8-13 Hz alpha peak prominence
f_alpha_idx = freq >= 8 & freq <= 13;
f_flank_idx = (freq >= 6 & freq < 8) | (freq > 13 & freq <= 15);

alpha_peak_db  = max(psd_db(cand_ecg_ics, f_alpha_idx), [], 2);
baseline_floor = mean(psd_db(cand_ecg_ics, f_flank_idx), 2);

% If alpha power exceeds the surrounding floor by > 3 dB, it contains brain rhythm
has_alpha_peak = (alpha_peak_db - baseline_floor) > 6.0;

% Quantify harmonic comb ripple in the 2-12 Hz window
f_sub = freq >= 1 & freq <= 12;
p_sub = psd_db(cand_ecg_ics, f_sub);

max_freq_ecg = 35;
has_ecg_spectrum = false(size(has_alpha_peak));
for i_comp = 1:n_cands
    ic_idx = cand_ecg_ics(i_comp);

    % Detrend gross slope across 2-12 Hz to isolate ripples
    % p_ripple = detrend(p_sub(i_comp, :));
    p_ripple = p_sub(i_comp, :);

    % Detect harmonic teeth
    num_peaks = length(findpeaks(p_ripple, 'MinPeakProminence', 2));
    has_cardiac_comb = (num_peaks > 2);

    % High-frequency power attenuation check (HF/LF < 10%)
    hf_lf_ratio = mean(10.^(psd_db(ic_idx, freq >= (max_freq_ecg + 5)) / 10)) / ...
        mean(10.^(psd_db(ic_idx, freq >= 1 & freq <= max_freq_ecg) / 10)) * 100;
    is_low_hf = (hf_lf_ratio < 35);

    % Combined spectral check for genuine ECG
    has_ecg_spectrum(i_comp) = has_cardiac_comb && is_low_hf;
end

% Save for checks
if any(has_ecg_spectrum)
    psd_db_ecg = psd_db(has_ecg_spectrum, :);
    path_save  = 'C:\DATA\MATLAB\myCodes\preprocessing\files\weights\ecg';
    if ~exist(path_save, 'dir'), mkdir(path_save); end
    file_name  = [EEG.ALSUTRECHT.subject.id '_heart_spectra.mat'];
    save(fullfile(path_save, file_name), 'freq', 'psd_db_ecg');
end

% -------------------------------------------------------------------------
% 4. Visualisation Setup (2-Row Nested Tiled Layout)
% -------------------------------------------------------------------------
if cfg.figure.plot
    if n_cands <= 8
        colours = brewermap(max(3, n_cands), 'Dark2');
    elseif n_cands <= 12
        colours = brewermap(n_cands, 'Paired');
    else
        colours = lines(n_cands);
    end

    fig_width  = max(950, 260 * n_cands + 100);
    fig_height = 680;
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.figure.visible);

    % Master layout spanning 100% of figure area
    t_main = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

    % Row 1: 1 row x 3 columns (Cols 1-2: Power Spectra; Col 3: Dedicated Legend)
    t_top = tiledlayout(t_main, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_top.Layout.Tile = 1;
    ax_psd = nexttile(t_top, [1, 2]);
    hold(ax_psd, 'on');

    % Shaded region highlighting typical EMG band (> 30 Hz)
    y_bounds = [min(psd_db(:)) - 3, max(psd_db(:)) + 3];
    patch(ax_psd, [max_freq_ecg, cfg.max_freq, cfg.max_freq, max_freq_ecg], ...
        [y_bounds(1), y_bounds(1), y_bounds(2), y_bounds(2)], ...
        [0.94, 0.94, 0.94], 'EdgeColor', 'none', 'DisplayName', 'EMG Band');

    % Boundary marker
    xline(ax_psd, max_freq_ecg, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Background: All components in faint grey
    h_all = plot(ax_psd, freq, psd_db', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);

    % Candidate components
    h_cands = gobjects(n_cands, 1);
    leg_txt = cell(n_cands, 1);

    for i = 1:n_cands
        ic_idx = cand_ecg_ics(i);
        h_cands(i) = plot(ax_psd, freq, psd_db(ic_idx, :), ...
            'Color', colours(i, :), 'LineWidth', 2.0);

        if has_alpha_peak(i)
            alpha_str = 'Yes';
        else
            alpha_str = 'No';
        end

        if has_ecg_spectrum(i)
            ecg_str = 'Yes';
        else
            ecg_str = 'No';
        end

        leg_txt{i} = sprintf('IC%2d (Alpha: %3s | ECG: %3s)', ic_idx, alpha_str, ecg_str);
    end

    xlim(ax_psd, [1, cfg.max_freq]);
    ylim(ax_psd, y_bounds);
    grid(ax_psd, 'on');
    set(ax_psd, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
    xlabel(ax_psd, 'Frequency (Hz)', 'FontWeight', 'bold');
    ylabel(ax_psd, 'Power (10\cdotlog_{10} \muV^2/Hz)', 'FontWeight', 'bold');
    title(ax_psd, 'Power Spectra of Candidate ECG Components', 'FontSize', 12, 'FontWeight', 'bold');
    subtitle(ax_psd, 'ECG has harmonic comb ripple (<15 Hz) & belly spectrum (<30 Hz)', ...
        'FontAngle', 'italic', 'Color', [0.4, 0.4, 0.4]);

    % Dock legend inside Tile 3 of the top sub-layout
    lgd = legend(ax_psd, [h_all(1); h_cands], [{'All ICs'}; leg_txt], ...
        'Box', 'off', 'FontSize', 8);
    lgd.Layout.Tile = 3;
    hold(ax_psd, 'off');

    % ---------------------------------------------------------------------
    % Row 2: Topoplots of Candidate ICs (Full Bottom Width)
    % ---------------------------------------------------------------------
    t_bot = tiledlayout(t_main, 1, n_cands, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_bot.Layout.Tile = 2;

    winv = EEG.icawinv;
    for i_comp = 1:n_cands
        ax_topo = nexttile(t_bot);
        ic_idx     = cand_ecg_ics(i_comp);
        pc_weights = winv(:, ic_idx);
        cmax       = max(abs(pc_weights));
        if cmax == 0, cmax = 1; end

        mytopoplot(pc_weights, [], '', ax_topo, [-cmax, cmax]);

        % Highlight when component lacks alpha and possesses genuine ECG spectral properties
        is_ecg = ~has_alpha_peak(i_comp) && has_ecg_spectrum(i_comp);

        if is_ecg
            title_str = sprintf('\\bf★ IC%d (ECG) ★', ic_idx);
            % title_col = [0.85, 0.15, 0.15]; % Signal red highlight
            title_col = colours(i_comp, :);

            set(ax_topo, 'Box', 'on', ...
                'XColor', title_col, ...
                'YColor', title_col, ...
                'LineWidth', 2.5);
        else
            title_str = sprintf('IC%d', ic_idx);
            title_col = colours(i_comp, :);
            set(ax_topo, 'Box', 'off');
        end

        title(ax_topo, title_str, 'Color', title_col, 'FontWeight', 'bold', 'FontSize', 11);
    end

    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_ecg_spectra'], [22, 14]);
end
end
function bad_ic = is_false_heog(bad_ic, ICs, templates_ica, ICLabel)
if isempty(bad_ic); return; end

blink_templates = [templates_ica.Blinkweights0, templates_ica.Blinkweights1];
max_sim_blink = calc_max_sim(ICs(:, bad_ic), blink_templates, 'cosine');

false_is_blink = max_sim_blink > 0.85;
false_is_brain = is_likely_brain(bad_ic, ICLabel, 0.60);

false_all = false_is_blink(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[HEOG Validation] Pruned %d/%d false-positive HEOG ICs.\n', sum(false_all), length(false_all));
end

function bad_ic = is_false_emg(bad_ic, ICs, templates_ica, ICLabel)
if isempty(bad_ic); return; end

blink_fields = fieldnames(templates_ica);
blink_mask   = contains(lower(blink_fields), 'blink');
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

heog_mask = contains(lower(blink_fields), 'saccade');
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

false_is_brain = is_likely_brain(bad_ic, ICLabel, 0.60);
false_all = false_is_blink(:) | false_is_heog(:) | false_is_brain(:);
bad_ic(false_all) = [];

fprintf('[EMG Validation]  Pruned %d/%d false-positive EMG ICs (Blink: %d, HEOG: %d, Brain: %d).\n', ...
    sum(false_all), length(false_all), sum(false_is_blink), sum(false_is_heog), sum(false_is_brain));
end

function is_brain = is_likely_brain(bad_ic, ICLabel, prob_threshold)
is_brain = (ICLabel.cvec(bad_ic) == 1) & (ICLabel.pvec(bad_ic) >= prob_threshold);
end

function [bad_ics, sim_scores, dist_scores, report] = match_ica_template(icawinv, chanlocs, templates_ica, template_type, cfg)
if nargin < 5; cfg = struct(); end

assert(~isempty(icawinv) && ismatrix(icawinv), 'icawinv must be a non-empty 2D matrix.');
icawinv = double(icawinv);
[num_chans, num_ics] = size(icawinv);

if ~isempty(chanlocs) && isstruct(chanlocs)
    chan_labels = upper({chanlocs.labels});
else
    chan_labels = {};
end

all_field_names = fieldnames(templates_ica);
matched_idx = contains(lower(all_field_names), lower(template_type));
if ~any(matched_idx) && strcmpi(template_type, 'saccade')
    matched_idx = contains(lower(all_field_names), 'horizontal');
end

assert(any(matched_idx), 'No fields matching "%s" found in templates_ica struct.', template_type);
matched_fields = all_field_names(matched_idx);

W_tgt = [];
for i_f = 1:length(matched_fields)
    current_template = double(templates_ica.(matched_fields{i_f}));
    assert(size(current_template, 1) == num_chans, ...
        'Field %s has %d channels; expected %d.', matched_fields{i_f}, size(current_template, 1), num_chans);
    W_tgt = [W_tgt, current_template]; %#ok<AGROW>
end

num_variants = size(W_tgt, 2);

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

[sim_scores, sim_matrix, best_variant_idx] = calc_max_sim(icawinv, W_tgt, cfg.metric);

dist_scores = NaN(1, num_ics);
pass_dist   = true(1, num_ics);

if cfg.check_distribution
    total_energy = sum(icawinv.^2, 1) + eps;
    switch lower(template_type)
        case 'blink'
            roi_mask = ismember(chan_labels, cfg.biosemi_blink);
            if any(roi_mask)
                dist_scores = sum(icawinv(roi_mask, :).^2, 1) ./ total_energy;
                pass_dist   = dist_scores >= cfg.dist_thresh;
            end
        case 'saccade'
            roi_mask = ismember(chan_labels, cfg.biosemi_saccade);
            if any(roi_mask)
                dist_scores = sum(icawinv(roi_mask, :).^2, 1) ./ total_energy;
                pass_dist   = dist_scores >= cfg.dist_thresh;
            end
        case 'heart'
            dist_scores = max(icawinv.^2, [], 1) ./ total_energy;
            pass_dist   = dist_scores <= cfg.dist_thresh;
    end
else
    cfg.dist_thresh = [];
end

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
end

function [max_sim, sim_matrix, best_variant_idx] = calc_max_sim(W_ic, W_tgt, metric)
if nargin < 3 || isempty(metric); metric = 'cosine'; end

W_ic  = double(W_ic);
W_tgt = double(W_tgt);

switch lower(metric)
    case 'cosine'
        norm_ic   = sqrt(sum(W_ic.^2, 1));
        norm_tgt  = sqrt(sum(W_tgt.^2, 1))';
        dot_prods = abs(W_tgt' * W_ic);
        sim_matrix = dot_prods ./ (norm_tgt * norm_ic + eps);

    case 'pearson'
        sim_matrix = abs(corr(W_ic, W_tgt, 'type', 'Pearson'))';

    otherwise
        error('Unknown metric "%s". Use ''cosine'' or ''pearson''.', metric);
end

[max_sim, best_variant_idx] = max(sim_matrix, [], 1);
end