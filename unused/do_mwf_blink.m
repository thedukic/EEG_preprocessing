function [EEG, W, mask_blink] = do_mwf_blink(EEG, cfg)
% DO_MWF_BLINK Detects ocular excursions on continuous 2D EEG using the dedicated
% VEOG channel, attenuates blink artefacts on EEG channels via spatio-temporal MWF,
% and renders topographic verification of power changes using mytopoplot.

if nargin < 2, cfg = struct(); end

% 1. Validation & Default Configuration
assert(ismatrix(EEG.data), 'EEG.data must be 2D continuous data. Run this function prior to epoching.');

if ~isfield(cfg, 'veog_chan'),        cfg.veog_chan        = 'VEOG'; end
if ~isfield(cfg, 'iqr_mult'),         cfg.iqr_mult         = 3.0;    end
if ~isfield(cfg, 'min_deflect'),      cfg.min_deflect      = 50;     end
if ~isfield(cfg, 'pad_sec'),          cfg.pad_sec          = 0.00;   end
% if ~isfield(cfg, 'delay_ms'),         cfg.delay_ms         = 10;     end
% if ~isfield(cfg, 'delay_spacing_ms'), cfg.delay_spacing_ms = 16;     end
if ~isfield(cfg, 'do_plot'),          cfg.do_plot          = true;   end

fs = EEG.srate;

% 2. VEOG Channel Extraction
% Locate VEOG channel
if ischar(cfg.veog_chan) || isstring(cfg.veog_chan)
    idx_veog = find(strcmpi({EEG.chanlocs.labels}, cfg.veog_chan), 1);
else
    idx_veog = cfg.veog_chan;
end
assert(~isempty(idx_veog), 'VEOG channel [%s] not found in DATA.', string(cfg.veog_chan));

% Extract VEOG trace and remove slow baseline drift
eye_sig = double(EEG.data(idx_veog, :));
eye_sig = eye_sig - movmedian(eye_sig, fs);

% 3. Non-Parametric Threshold Detection
veog_IQR   = iqr(eye_sig);
veog_75P   = prctile(eye_sig, 75);

det_thresh = max(cfg.min_deflect, veog_75P + (cfg.iqr_mult * veog_IQR));
raw_hits   = eye_sig > det_thresh;

% Dilate mask symmetrically around each detected point
pad_samples = round(cfg.pad_sec * fs);
kernel      = ones(1, 2 * pad_samples + 1);
mask_blink  = conv(double(raw_hits), kernel, 'same') > 0;

n_art_samples   = sum(mask_blink);
n_clean_samples = sum(~mask_blink);

fprintf('--- Continuous MWF Blink Cleaning ---\n');
fprintf('  Detector Channel: %s\n', cfg.veog_chan);
fprintf('  Threshold:        %.2f uV\n', det_thresh);
fprintf('  Artefact masked:  %d samples (%.2f s)\n', n_art_samples, n_art_samples / fs);
fprintf('  Clean retained:   %d samples (%.2f s)\n', n_clean_samples, n_clean_samples / fs);

% Initialise logging struct
EEG.ALSUTRECHT.mwf.status                 = false;
EEG.ALSUTRECHT.mwf.signalToErrorRatio     = NaN;
EEG.ALSUTRECHT.mwf.artifactToResidueRatio = NaN;
W = [];

% Minimum data checks
min_req_samples = round(2.0 * fs);
if n_art_samples < min_req_samples
    warning('Insufficient blink data detected (< 2.0 s). Skipping MWF.');
    return;
end

if n_clean_samples < min_req_samples
    warning('Insufficient clean baseline segments remaining. Skipping MWF.');
    return;
end

% 4. MWF Configuration & Dimension Reduction
mask_eeg = strcmpi({EEG.chanlocs.type}, 'EEG');
raw_eeg  = double(EEG.data(mask_eeg, :));

% Sample 30s to determine true rank (accounts for avg ref & interpolated chans)
data_rank = rank(raw_eeg);

% Project down to full-rank orthogonal subspace via SVD/PCA
[U, ~, ~] = svd(raw_eeg, 'econ');
proj_mat  = U(:, 1:data_rank)';

% Transform continuous data to PCA component space
comp_data = proj_mat * raw_eeg;

% Configure and run MWF on the compressed subspace
% params = mwf_params('delay', delay_samples, 'delay_spacing', spacing_samples);

% 1. Pure Spatial Filter (Zero temporal lag for volume-conducted blinks)
% 2. Retain strictly the 1 or 2 dominant generalized eigenvalues
% params = mwf_params(...
%     'delay',   0, ...
%     'rank',    'first', ...
%     'rankopt', 1); % 1 captures the vertical dipole; 2 catches Bell's phenomenon / lateral tilt

params = mwf_params(...
    'delay', 10, ...
    'delay_spacing', 5, ...
    'rank',    'first', ...
    'rankopt', 15);

% Run MWF directly on the full-rank PCA subspace or raw EEG
lastwarn('');
[comp_clean, ~, W, SER, ARR] = mwf_process(comp_data, mask_blink, params);

% Project back to original scalp sensor space
eeg_clean = proj_mat' * comp_clean;

% % Visual inspection
% EEG_NEW = EEG;
% EEG_NEW.data(mask_eeg, :) = eeg_clean;
% vis_artifacts(EEG_NEW, EEG);

% 5. Topographic Power Diagnostics
if cfg.do_plot
    % Demean prior to power computation to isolate variance
    raw_art   = raw_eeg(:, mask_blink);
    clean_art = eeg_clean(:, mask_blink);

    raw_art   = raw_art - mean(raw_art, 2);
    clean_art = clean_art - mean(clean_art, 2);

    % Power during blink segments (mean squared signal)
    pow_pre  = mean(raw_art.^2, 2);
    pow_post = mean(clean_art.^2, 2);

    % Attenuation in decibels
    delta_db_blink = 10 * log10(pow_post ./ pow_pre);

    % Baseline preservation check during clean segments
    raw_clean_seg   = raw_eeg(:, ~mask_blink);
    clean_clean_seg = eeg_clean(:, ~mask_blink);
    raw_clean_seg   = raw_clean_seg - mean(raw_clean_seg, 2);
    clean_clean_seg = clean_clean_seg - mean(clean_clean_seg, 2);

    pow_clean_pre  = mean(raw_clean_seg.^2, 2);
    pow_clean_post = mean(clean_clean_seg.^2, 2);
    delta_db_clean = 10 * log10(pow_clean_post ./ pow_clean_pre);

    % Render 1x4 diagnostic summary
    fh = figure('Color', 'w', 'Position', [100, 100, 1400, 340], 'Name', 'MWF Blink Power Verification');
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    max_raw_pow = max(pow_pre);
    % min_db      = floor(min(delta_db_blink(isfinite(delta_db_blink))));
    % if isempty(min_db), min_db = -10; end

    % 1. Pre-Cleaning Blink Power
    ax1 = nexttile;
    mytopoplot(pow_pre, [], '', ax1, [0, max_raw_pow]);
    hcb1 = colorbar(ax1);
    hcb1.Title.String = '\muV^2';
    title(ax1, sprintf('Pre-MWF Blink Power\n(Max: %.1f \\muV^2)', max_raw_pow), 'FontWeight', 'bold');

    % 2. Post-Cleaning Blink Power (Matched Color Limits)
    ax2 = nexttile;
    mytopoplot(pow_post, [], '', ax2, [0, max_raw_pow]);
    hcb2 = colorbar(ax2);
    hcb2.Title.String = '\muV^2';
    title(ax2, sprintf('Post-MWF Blink Power\n(Max: %.1f \\muV^2)', max(pow_post)), 'FontWeight', 'bold');

    % 3. Blink Attenuation (dB)
    ax3 = nexttile;
    mytopoplot(delta_db_blink, [], '', ax3);
    hcb3 = colorbar(ax3);
    hcb3.Title.String = 'dB';
    title(ax3, sprintf('Blink Attenuation\n(Min: %.1f dB)', min(delta_db_blink)), 'FontWeight', 'bold');

    % 4. Clean Baseline Distortion Check (dB)
    ax4 = nexttile;
    mytopoplot(delta_db_clean, [], '', ax4);
    hcb4 = colorbar(ax4);
    hcb4.Title.String = 'dB';
    title(ax4, sprintf('Clean Preservation\n(Mean: %.2f dB)', mean(delta_db_clean)), 'FontWeight', 'bold');

    % Save figure if subject directory exists
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_mwf_blink_power'], [36 10]);

end

% 6. Validate Output
flag_successful = true;
[warn_msg, ~]   = lastwarn;

if contains(warn_msg, "eigenvectors") || contains(warn_msg, "singular") || isnan(SER) || isnan(ARR) || ~isreal(eeg_clean)
    warning('Something went wrong with MWF: %s', warn_msg);
    flag_successful = false;
end

% 7. Store Outputs and Metrics
if flag_successful
    EEG.data(mask_eeg, :) = eeg_clean;
    EEG.ALSUTRECHT.mwf.status                 = true;
    EEG.ALSUTRECHT.mwf.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.mwf.artifactToResidueRatio = ARR;

    fprintf('Signal-to-error ratio:     %1.2f\n', SER);
    fprintf('Artifact-to-residue ratio: %1.2f\n', ARR);
end

end