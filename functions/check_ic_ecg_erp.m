function [is_ecg, stats, fh] = check_ic_ecg_erp(EEG, EXT, ECGepochs, cfg)
% CHECK_IC_ECG_ERP
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

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 4, cfg = struct(); end

if ~isfield(cfg, 'epoch_window_ms'), cfg.epoch_window_ms = [-200, 250]; end
if ~isfield(cfg, 'min_corr'),        cfg.min_corr        = 0.80;        end
if ~isfield(cfg, 'min_snr'),         cfg.min_snr         = 3.5;         end
if ~isfield(cfg, 'max_lag_ms'),      cfg.max_lag_ms      = 30;          end
if ~isfield(cfg, 'do_plot'),         cfg.do_plot         = true;        end

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
else
    ch_idx = 1:size(EEG.icaweights, 2);
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
else
    ecg_idx = find(strcmpi({EXT.chanlocs.labels}, 'ECG')  | ...
        strcmpi({EXT.chanlocs.labels}, 'EKG')  | ...
        strcmpi({EXT.chanlocs.labels}, 'EXG1') | ...
        strcmpi({EXT.chanlocs.labels}, 'BIP1') | ...
        strcmpi({EXT.chanlocs.labels}, 'HR'), 1);
    if isempty(ecg_idx)
        warning('ECG label not found in EXT. Defaulting to channel 1 of EXT.');
        ecg_idx = 1;
    end
end

n_pnts_ext    = EXT.pnts * EXT.trials;
raw_ecg       = double(reshape(EXT.data(ecg_idx, :, :), 1, n_pnts_ext));
total_samples = min(n_pnts_eeg, n_pnts_ext);

ica_act = ica_act(:, 1:total_samples);
raw_ecg = raw_ecg(1:total_samples);

% -------------------------------------------------------------------------
% 3. Extract Epoched Data and Compute Trial-Averaged ERPs
% -------------------------------------------------------------------------
n_epochs   = size(ECGepochs, 1);
epoch_len  = ECGepochs(1, 2) - ECGepochs(1, 1) + 1;
max_lag_sm = round((cfg.max_lag_ms / 1000) * fs);

ecg_trials = zeros(epoch_len, n_epochs);
ic_trials  = zeros(n_ics, epoch_len, n_epochs);

valid_epochs = 0;
for ep = 1:n_epochs
    i_start = ECGepochs(ep, 1);
    i_end   = ECGepochs(ep, 2);

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

% Baseline definition: First 80 ms of epoch (well before the QRS complex)
base_samples = max(2, round(0.080 * fs));
base_pts     = 1:base_samples;

ecg_erp = ecg_erp - mean(ecg_erp(base_pts));
for k = 1:n_ics
    ic_erps(k, :) = ic_erps(k, :) - mean(ic_erps(k, base_pts));
end

% -------------------------------------------------------------------------
% 4. Morphological Cross-Correlation & SNR Evaluation
% -------------------------------------------------------------------------
max_r    = zeros(n_ics, 1);
best_lag = zeros(n_ics, 1);
snr_ic   = zeros(n_ics, 1);
is_ecg   = false(n_ics, 1);

ecg_norm = (ecg_erp - mean(ecg_erp)) / (std(ecg_erp) * sqrt(epoch_len - 1));

for k = 1:n_ics
    cur_ic_erp = ic_erps(k, :)';

    % A. Signal-to-Noise Ratio (Peak-to-Peak vs Baseline SD)
    base_sd = std(cur_ic_erp(base_pts));
    if base_sd == 0, base_sd = eps; end
    snr_ic(k) = (max(cur_ic_erp) - min(cur_ic_erp)) / base_sd;

    % B. Cross-Correlation across physiological lag window
    ic_norm = (cur_ic_erp - mean(cur_ic_erp)) / (std(cur_ic_erp) * sqrt(epoch_len - 1));
    [r_lags, lags] = xcorr(ic_norm, ecg_norm, max_lag_sm, 'coeff');

    [peak_corr, max_idx] = max(abs(r_lags));
    max_r(k)    = peak_corr;
    best_lag(k) = (lags(max_idx) / fs) * 1000;

    % C. Strict Decision Rule
    if (max_r(k) >= cfg.min_corr) && (snr_ic(k) >= cfg.min_snr) && (abs(best_lag(k)) <= cfg.max_lag_ms)
        is_ecg(k) = true;
    end
end

% Optional Dominance Filter: If multiple ICs pass, keep only those close to the max correlation
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
    k = stats.flagged_ics(i);
    fprintf('   -> IC %2d: |r| = %.2f | Lag = %+5.1f ms | SNR = %4.1f\n', ...
        k, max_r(k), best_lag(k), snr_ic(k));
end
fprintf('==================================================\n\n');

% -------------------------------------------------------------------------
% 5. Visualisation
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot
    return;
end

time_vec = linspace(epoch_start_ms, epoch_end_ms, epoch_len);

flagged_ics = stats.flagged_ics;
n_flagged   = length(flagged_ics);
n_cols      = max(3, n_flagged);

fh = figure('Color', 'w');

% Top Panel: R-Peak Triggered ERP Traces
subplot(2, n_cols, 1:n_cols);
hold on;

% Background: All non-cardiac ICs (normalised)
ic_erps_norm = stats.ic_erps ./ max(abs(stats.ic_erps), [], 2);
h_all = plot(time_vec, ic_erps_norm', 'Color', [0.82, 0.82, 0.82], 'LineWidth', 0.8);

% Reference ECG ERP (from EXT)
h_ecg = plot(time_vec, stats.ecg_erp / max(abs(stats.ecg_erp)), 'k', 'LineWidth', 2.5);

% Flagged Cardiac ICs (polarity-aligned)
colours   = lines(max(7, n_flagged));
h_flagged = gobjects(n_flagged, 1);
leg_entries = {'All ICs (Background)', sprintf('Reference ECG (EXT.%s)', stats.ecg_channel)};

for i = 1:n_flagged
    ic_idx  = flagged_ics(i);
    ic_wave = stats.ic_erps(ic_idx, :);

    if corr(ic_wave', stats.ecg_erp) < 0
        ic_wave = -ic_wave;
    end
    ic_wave_norm = ic_wave / max(abs(ic_wave));

    h_flagged(i) = plot(time_vec, ic_wave_norm, 'Color', colours(i, :), 'LineWidth', 2.0);
    leg_entries{end+1} = sprintf('IC %d (|r|=%.2f, SNR=%.1f, lag=%+.1fms)', ...
        ic_idx, stats.max_r(ic_idx), stats.snr(ic_idx), stats.best_lag_ms(ic_idx));
end

xline(0, '--r', 'R-Peak (t = 0)', 'LineWidth', 1.2, 'LabelOrientation', 'aligned');
grid on;
xlim([epoch_start_ms, epoch_end_ms]);
ylim([-1.2, 1.2]);
xlabel('Time relative to R-peak (ms)');
ylabel('Normalised Potential (a.u.)');
title(sprintf('R-Peak Triggered QRS Average (%d Valid Epochs)', stats.valid_epochs), 'FontSize', 11);

if n_flagged > 0
    legend([h_all(1), h_ecg, h_flagged'], leg_entries, 'Location', 'northeast', 'FontSize', 8);
else
    legend([h_all(1), h_ecg], {'All ICs (Background)', sprintf('Reference ECG (EXT.%s - No ICs Flagged)', stats.ecg_channel)}, 'Location', 'northeast');
end
pbaspect([2.8 1 1]);

% Bottom Panel: Topoplots of Flagged Cardiac ICs
if n_flagged == 0
    sbh = subplot(2, n_cols, (n_cols + 1):(2 * n_cols));
    text(sbh, 0.5, 0.5, 'No Independent Components reached cardiac classification thresholds', ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontAngle', 'italic', 'Color', [0.4 0.4 0.4]);
    axis(sbh, 'off');
else
    for i = 1:min(n_cols, n_flagged)
        sbh = subplot(2, n_cols, n_cols + i);

        ic_idx     = flagged_ics(i);
        pc_weights = winv(:, ic_idx);
        cmax       = max(abs(pc_weights));
        if cmax == 0, cmax = 1; end

        title_str = sprintf('IC %d Topo (|r|=%.2f)', ic_idx, stats.max_r(ic_idx));
        mytopoplot(pc_weights, [], title_str, sbh, [-cmax, cmax]);
    end
end

drawnow;

end