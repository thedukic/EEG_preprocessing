function [is_ecg, stats, fh] = check_ecg_erp(EEG, EXT, ECGepochs, cfg)
% CHECK_ECG_ERP
% Identifies primary cardiac (ECG) ICs by comparing their QRS ERP morphology,
% SNR, and latency against the continuous reference ECG channel in EXT.
%
% Inputs:
%   EEG              : EEGLAB dataset structure (scalp EEG + ICA decomposition)
%   EXT              : EEGLAB dataset structure (external bipolar leads)
%   ECGepochs        : [N_trials x 2] Sample index matrix [start, end] per R-peak
%   cfg              : Configuration structure:
%                      cfg.epoch_window_ms : [start_ms, end_ms] relative to R-peak (default: [-200, 250])
%                      cfg.min_corr        : Min absolute correlation (default: 0.80)
%                      cfg.min_snr         : Min peak-to-peak to baseline SNR (default: 3.5)
%                      cfg.max_lag_ms      : Max allowed conduction delay in ms (default: 30)
%                      cfg.ecg_chan        : ECG channel label or index in EXT (default: auto-detect)
%                      cfg.do_plot         : Boolean to render diagnostic figure (default: true)
%
% Outputs:
%   is_ecg           : Logical column vector [n_ics x 1] (true for primary cardiac ICs)
%   stats            : Diagnostic statistics structure
%   fh               : Figure handle (or empty if cfg.do_plot = false)

% Expected Lag of a True Heart IC
% The expected lag depends entirely on the physical source of the cardiac artefact:
% - Electrical ECG Artefact: The propagation of the cardiac electrical field through the body to the scalp is virtually instantaneous.
%   The expected lag for these components is 0 to 10 ms.
% - Ballistocardiogram (BCG) Artefact: This represents the mechanical movement of the scalp and electrodes driven by the systolic blood pulse in the scalp arteries.
%   Because this relies on physical pulse transit time, the expected lag is 100 to 250 ms after the R-peak.

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 4, cfg = struct(); end

if ~isfield(cfg, 'epoch_window_ms'), cfg.epoch_window_ms = [-200, 250]; end
if ~isfield(cfg, 'min_corr'),        cfg.min_corr        = 0.80;        end
if ~isfield(cfg, 'min_snr'),         cfg.min_snr         = 3.5;         end
if ~isfield(cfg, 'max_lag_ms'),      cfg.max_lag_ms      = 30;          end
if ~isfield(cfg, 'min_ipr'),         cfg.min_ipr         = 12.0;        end
if ~isfield(cfg, 'max_top3'),        cfg.max_top3        = 0.45;        end
if ~isfield(cfg, 'min_extent'),      cfg.min_extent      = 0.08;        end
if ~isfield(cfg, 'do_plot'),         cfg.do_plot         = true;        end
if ~isfield(cfg, 'plot_visible'),    cfg.plot_visible    = 'on';        end

assert(EEG.srate == EXT.srate, 'Sampling rates of EEG (%d Hz) and EXT (%d Hz) must match.', EEG.srate, EXT.srate);
fs = EEG.srate;

% Compute exact sample offset for t = 0 (R-peak) based on epoch window
epoch_start_ms   = cfg.epoch_window_ms(1);
epoch_end_ms     = cfg.epoch_window_ms(2);
r_peak_samp_idx  = round(abs(epoch_start_ms) / 1000 * fs) + 1;

% -------------------------------------------------------------------------
% 2. Extract Data: ICA Activations (from EEG) & Reference ECG (from EXT)
% -------------------------------------------------------------------------
if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'ica') && isfield(EEG.ALSUTRECHT.ica, 'icachansind')
    ch_idx = EEG.ALSUTRECHT.ica.icachansind;
    w_mat  = EEG.ALSUTRECHT.ica.icaweights * EEG.ALSUTRECHT.ica.icasphere;
    winv   = EEG.ALSUTRECHT.ica.icawinv;
elseif isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    ch_idx = EEG.icachansind;
    w_mat  = EEG.icaweights * EEG.icasphere;
    winv   = EEG.icawinv;
end

n_pnts_eeg = EEG.pnts * EEG.trials;
scalp_data = double(reshape(EEG.data(ch_idx, :, :), length(ch_idx), n_pnts_eeg));
ica_act    = w_mat * scalp_data;
n_ics      = size(ica_act, 1);

% Locate reference ECG channel in EXT
if isfield(cfg, 'ecg_chan') && ~isempty(cfg.ecg_chan)
    if ischar(cfg.ecg_chan) || isstring(cfg.ecg_chan)
        ecg_idx = find(strcmpi({EXT.chanlocs.labels}, cfg.ecg_chan), 1);
    else
        ecg_idx = cfg.ecg_chan;
    end
end

n_pnts_ext    = EXT.pnts * EXT.trials;
raw_ecg       = double(reshape(EXT.data(ecg_idx, :, :), 1, n_pnts_ext));
total_samples = min(n_pnts_eeg, n_pnts_ext);

ica_act = ica_act(:, 1:total_samples);
raw_ecg = raw_ecg(1:total_samples);

% Apply a 1-25 Hz bandpass filter to both signals for morphological matching
[b, a] = butter(2, [1, 25]/(fs/2), 'bandpass');
ica_act = filtfilt(b, a, ica_act')';
raw_ecg = filtfilt(b, a, raw_ecg);

% -------------------------------------------------------------------------
% 3. Extract Epoched Data and Compute Trial-Averaged ERPs
% -------------------------------------------------------------------------
n_epochs   = size(ECGepochs, 1);
epoch_len  = ECGepochs(1, 2) - ECGepochs(1, 1) + 1;
max_lag_sm = round((cfg.max_lag_ms / 1000) * fs);

ecg_trials = zeros(epoch_len, n_epochs);
ic_trials  = zeros(n_ics, epoch_len, n_epochs);

valid_epochs = 0;
for i_epoch = 1:n_epochs
    i_start = ECGepochs(i_epoch, 1);
    i_end   = ECGepochs(i_epoch, 2);

    if i_start >= 1 && i_end <= total_samples
        valid_epochs = valid_epochs + 1;
        ecg_trials(:, valid_epochs)   = raw_ecg(i_start:i_end)';
        ic_trials(:, :, valid_epochs) = ica_act(:, i_start:i_end);
    end
end

ecg_trials = ecg_trials(:, 1:valid_epochs);
ic_trials  = ic_trials(:, :, 1:valid_epochs);

% Compute R-peak triggered averages
ecg_erp = mean(ecg_trials, 2);
ic_erps = mean(ic_trials, 3);

% Baseline definition: First 100 or 80 ms of epoch (well before the QRS complex)
base_samples = max(2, round(0.100 * fs));
base_pts     = 1:base_samples;

ecg_erp = ecg_erp - mean(ecg_erp(base_pts));
for i_comp = 1:n_ics
    ic_erps(i_comp, :) = ic_erps(i_comp, :) - mean(ic_erps(i_comp, base_pts));
end

% -------------------------------------------------------------------------
% Spatial Focality & Far-Field Validation
% -------------------------------------------------------------------------
winv     = EEG.icawinv;
abs_winv = abs(winv);

% 1. Inverse Participation Ratio (Effective Channel Count)
ipr = (sum(winv.^2, 1).^2) ./ sum(winv.^4, 1);

% 2. Top-3 Channel Energy Fraction
sort_sq_winv = sort(winv.^2, 1, 'descend');
top3_ratio   = sum(sort_sq_winv(1:3, :), 1) ./ sum(sort_sq_winv, 1);

% 3. Half-Maximum Spatial Extent (Fraction of montage >= 50% of peak)
peak_vals  = max(abs_winv, [], 1);
extent_50  = mean(abs_winv >= (0.5 * peak_vals), 1);

% % 4. Rim Proximity Check (Optional but highly recommended)
% % Flags whether the maximum absolute weight lies on the outermost rim
% if isfield(EEG.chanlocs, 'radius')
%     % EEGLAB polar radius: rim channels typically have radius >= 0.48
%     is_rim_chan = [EEG.chanlocs.radius] >= 0.48;
%     [~, max_ch] = max(abs_winv, [], 1);
%     peak_on_rim = is_rim_chan(max_ch);
% else
%     peak_on_rim = false(1, size(winv, 2));
% end

% -------------------------------------------------------------------------
% Decision Gate for True Far-Field Components (e.g. ECG)
% -------------------------------------------------------------------------
% True ECG must engage a broad spatial montage (IPR >= 12, Extent >= 8%)
% and cannot have > 45% of its power locked into just 3 channels.
is_far_field = (ipr >= cfg.min_ipr) & ...
    (top3_ratio <= cfg.max_top3) & ...
    (extent_50 >= cfg.min_extent);

% % If an IC is borderline but peaks right on the rim, reject it as focal EMG
% is_far_field = is_far_field & ~peak_on_rim;

% -------------------------------------------------------------------------
% 4. Morphological Cross-Correlation & SNR Evaluation
% -------------------------------------------------------------------------
max_r    = zeros(n_ics, 1);
best_lag = zeros(n_ics, 1);
snr_ic   = zeros(n_ics, 1);
is_ecg   = false(n_ics, 1);

ecg_norm = (ecg_erp - mean(ecg_erp)) / (std(ecg_erp) * sqrt(epoch_len - 1));

for i_comp = 1:n_ics
    cur_ic_erp = ic_erps(i_comp, :)';

    % A. Signal-to-Noise Ratio (Peak-to-Peak vs Baseline SD)
    base_sd = std(cur_ic_erp(base_pts));
    if base_sd == 0, base_sd = eps; end
    snr_ic(i_comp) = (max(cur_ic_erp) - min(cur_ic_erp)) / base_sd;

    % B. Cross-Correlation across physiological lag window
    ic_norm = (cur_ic_erp - mean(cur_ic_erp)) / (std(cur_ic_erp) * sqrt(epoch_len - 1));
    [r_lags, lags] = xcorr(ic_norm, ecg_norm, max_lag_sm, 'coeff');

    [peak_corr, max_idx] = max(abs(r_lags));
    max_r(i_comp)    = peak_corr;
    best_lag(i_comp) = (lags(max_idx) / fs) * 1000;

    % C. Strict Decision Rule
    if (max_r(i_comp) >= cfg.min_corr) && (snr_ic(i_comp) >= cfg.min_snr) && (abs(best_lag(i_comp)) <= cfg.max_lag_ms) && is_far_field(i_comp)
        is_ecg(i_comp) = true;
    end
end

% If multiple ICs pass, keep only those close to the max correlation
if sum(is_ecg) > 1
    top_r = max(max_r(is_ecg));
    % Reject components that are substantially weaker than the primary cardiac IC
    is_ecg(is_ecg & (max_r < (top_r - 0.08))) = false;
end

% Pack statistics structure
stats.ecg_erp          = ecg_erp;
stats.ic_erps          = ic_erps;
stats.max_r            = max_r;
stats.best_lag_ms      = best_lag;
stats.snr              = snr_ic;
stats.valid_epochs     = valid_epochs;
stats.flagged_ics      = find(is_ecg);
stats.ecg_channel      = EXT.chanlocs(ecg_idx).labels;
stats.r_peak_samp_idx  = r_peak_samp_idx;

% Console Summary
fprintf('\n==================================================\n');
fprintf('  QRS ERP Analysis: EXT [%s] vs IC Activations \n', stats.ecg_channel);
fprintf('==================================================\n');
fprintf('  Epoch Window: [%.1f, %.1f] ms | Analyzed Epochs: %d\n', epoch_start_ms, epoch_end_ms, valid_epochs);
fprintf('  Detection Criteria: |r| >= %.2f | SNR >= %.1f | Lag <= +/- %d ms\n', ...
    cfg.min_corr, cfg.min_snr, cfg.max_lag_ms);
fprintf('--------------------------------------------------\n');
fprintf('  Flagged Cardiac Components: %d\n', length(stats.flagged_ics));
for i = 1:length(stats.flagged_ics)
    i_comp = stats.flagged_ics(i);
    fprintf('   -> IC %2d: |r| = %.2f | Lag = %+5.1f ms | SNR = %4.1f\n', ...
        i_comp, max_r(i_comp), best_lag(i_comp), snr_ic(i_comp));
end
fprintf('==================================================\n\n');

% -------------------------------------------------------------------------
% 5. Visualisation
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot
    return;
end

time_vec    = linspace(epoch_start_ms, epoch_end_ms, epoch_len);
flagged_ics = stats.flagged_ics;
n_flagged   = length(flagged_ics);

% Distinct colour palette from brewermap Set2 (minimum 3 colours for ColorBrewer)
colours = brewermap(max(3, n_flagged), 'Set2');

% -------------------------------------------------------------------------
% Dynamic Layout Configuration
% -------------------------------------------------------------------------
if n_flagged == 0
    fh = figure('Color', 'w', 'Position', [100, 100, 950, 420], 'Visible', cfg.plot_visible);
    t_top = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_erp = nexttile(t_top, [1, 2]);
else
    fig_width  = max(950, 260 * n_flagged + 100);
    fig_height = 680;
    fh = figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height], 'Visible', cfg.plot_visible);

    % Master layout: 2 rows spanning 100% of figure width
    t_main = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

    % Row 1: 1 row x 3 columns (Cols 1-2: ERP traces; Col 3: Dedicated Legend)
    t_top = tiledlayout(t_main, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_top.Layout.Tile = 1;
    ax_erp = nexttile(t_top, [1, 2]);
end

% -------------------------------------------------------------------------
% Panel 1: R-Peak Triggered ERP Traces (Top Left)
% -------------------------------------------------------------------------
hold(ax_erp, 'on');

% Background: Non-flagged ICs (normalised)
ic_erps_norm = stats.ic_erps ./ max(abs(stats.ic_erps), [], 2);
h_all = plot(ax_erp, time_vec, ic_erps_norm', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);

% Reference ECG ERP (EXT lead)
h_ecg = plot(ax_erp, time_vec, stats.ecg_erp / max(abs(stats.ecg_erp)), 'k', 'LineWidth', 2.2);

% Flagged Cardiac ICs (polarity-aligned)
h_flagged   = gobjects(n_flagged, 1);
leg_entries = {'All ICs (Background)', sprintf('Reference ECG (EXT.%s)', stats.ecg_channel)};

for i_comp = 1:n_flagged
    ic_idx  = flagged_ics(i_comp);
    ic_wave = stats.ic_erps(ic_idx, :);

    if corr(ic_wave', stats.ecg_erp) < 0
        ic_wave = -ic_wave;
    end
    ic_wave_norm = ic_wave / max(abs(ic_wave));

    h_flagged(i_comp) = plot(ax_erp, time_vec, ic_wave_norm, 'Color', colours(i_comp, :), 'LineWidth', 2.0);
    leg_entries{end+1} = sprintf('IC%2d (|r| = %1.2f, SNR = %2.0f, lag = %+2.1f ms)', ...
        ic_idx, stats.max_r(ic_idx), stats.snr(ic_idx), stats.best_lag_ms(ic_idx));
end

grid(ax_erp, 'on');
set(ax_erp, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 10, 'GridAlpha', 0.4);
xlim(ax_erp, [epoch_start_ms, epoch_end_ms]);
ylim(ax_erp, [-1.3, 1.3]);
xlabel(ax_erp, 'Time relative to R-peak (ms)', 'FontWeight', 'bold');
ylabel(ax_erp, 'Normalised Potential (a.u.)', 'FontWeight', 'bold');
title(ax_erp, sprintf('R-Peak Triggered QRS Average (%d Valid Epochs)', stats.valid_epochs), ...
    'FontSize', 12, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Dedicated Legend Placement (Top Right Tile 3)
% -------------------------------------------------------------------------
if n_flagged == 0
    subtitle(ax_erp, 'No independent components reached cardiac classification thresholds', ...
        'FontAngle', 'italic', 'Color', [0.45, 0.45, 0.45], 'FontSize', 10);
    lgd = legend(ax_erp, [h_all(1), h_ecg], leg_entries, 'Box', 'off', 'FontSize', 8);
else
    % Omit 'Location' completely so it docks strictly inside Tile 3
    lgd = legend(ax_erp, [h_all(1), h_ecg, h_flagged'], leg_entries, ...
        'Box', 'off', 'FontSize', 8);
end

lgd.Layout.Tile = 3;
hold(ax_erp, 'off');

% -------------------------------------------------------------------------
% Panel 2: Topoplots of Flagged Cardiac ICs (Full Bottom Width)
% -------------------------------------------------------------------------
if n_flagged > 0
    % 1. Rank primarily by |R| with a flat tolerance dead-zone
    lag_tol   = 15; % ms: zero penalty within +/- 15 ms
    sigma_lag = 35; % ms: gentle rolloff beyond tolerance

    comp_scores = zeros(1, n_flagged);

    for i_comp = 1:n_flagged
        ic_idx  = flagged_ics(i_comp);
        r_val   = abs(stats.max_r(ic_idx));
        lag_val = abs(stats.best_lag_ms(ic_idx));

        % Flat-top penalty: 1.0 inside tolerance, soft exponential decay outside
        excess_lag = max(0, lag_val - lag_tol);
        lag_weight = exp(-(excess_lag^2) / (2 * sigma_lag^2));

        comp_scores(i_comp) = r_val * lag_weight;
    end

    [~, best_comp_pos] = max(comp_scores);

    % 2. Sub-layout spanning the full width of Row 2
    t_bot = tiledlayout(t_main, 1, n_flagged, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_bot.Layout.Tile = 2;

    for i_comp = 1:n_flagged
        ax_topo = nexttile(t_bot);
        ic_idx     = flagged_ics(i_comp);
        pc_weights = winv(:, ic_idx);
        cmax       = max(abs(pc_weights));
        if cmax == 0, cmax = 1; end

        mytopoplot(pc_weights, [], '', ax_topo, [-cmax, cmax]);

        % Highlight the primary cardiac component
        if i_comp == best_comp_pos
            title_str = sprintf('\\bf★ IC%d (BEST) ★', ic_idx);
            % title_col = [0.85, 0.15, 0.15];
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
end
end