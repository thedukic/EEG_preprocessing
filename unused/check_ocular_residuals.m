function [ocular_stats, fh] = check_ocular_residuals(DATA, cfg)
% CHECK_OCULAR_RESIDUALS
% Evaluates residual vertical (VEOG) and horizontal (HEOG) ocular artifact
% leakage in cleaned, epoched EEG data.
%
% Integrates:
%   - Epoch-level blink precedence (blinks suppress false-positive saccades)
%   - Physical-floor IQR event gating
%   - Delta-R contrast metrics and within-epoch ERP extraction
%   - 2-row diagnostic visualization (Row 1: Pre-Pruning, Row 2: Post-Pruning)
%
% Inputs:
%   DATA : Cleaned, epoched EEGLAB dataset (scalp channels + bipolar EOG)
%   cfg  : (Optional) Configuration structure:
%          cfg.veog_chan        : Label or index of VEOG in DATA (default: 'VEOG')
%          cfg.heog_chan        : Label or index of HEOG in DATA (default: 'HEOG')
%          cfg.max_corr         : Max acceptable residual correlation (default: 0.12)
%          cfg.trIQRsaccade     : IQR multiplier for HEOG saccades (default: 1.5)
%          cfg.trIQRblink       : IQR multiplier for VEOG blinks (default: 3.0)
%          cfg.min_phys_blink   : Minimum physical amplitude floor for blinks (default: 50 uV)
%          cfg.min_phys_saccade : Minimum physical step floor for saccades (default: 35 uV)
%          cfg.win_sec          : Half-window for ERP snippet extraction (default: 0.150 s)
%          cfg.do_plot          : Render diagnostic summary figure (default: true)
%
% Outputs:
%   ocular_stats : Struct containing pre/post r-maps, Delta-R maps, ERP metrics, and clean trial indices
%   fh           : Figure handle

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'veog_chan'),        cfg.veog_chan        = 'VEOG'; end
if ~isfield(cfg, 'heog_chan'),        cfg.heog_chan        = 'HEOG'; end
if ~isfield(cfg, 'max_corr'),         cfg.max_corr         = 0.12;   end
if ~isfield(cfg, 'trIQRsaccade'),     cfg.trIQRsaccade     = 1.5;    end
if ~isfield(cfg, 'trIQRblink'),       cfg.trIQRblink       = 3.0;    end
if ~isfield(cfg, 'min_phys_blink'),   cfg.min_phys_blink   = 50;     end % uV
if ~isfield(cfg, 'min_phys_saccade'), cfg.min_phys_saccade = 35;     end % uV
if ~isfield(cfg, 'win_sec'),          cfg.win_sec          = 0.150;  end % s
if ~isfield(cfg, 'do_plot'),          cfg.do_plot          = true;   end

fs = DATA.srate;

% -------------------------------------------------------------------------
% 2. Extract Scalp Data & Bipolar EOG Signals
% -------------------------------------------------------------------------
if isfield(DATA, 'ALSUTRECHT') && isfield(DATA.ALSUTRECHT, 'ica') && isfield(DATA.ALSUTRECHT.ica, 'icachansind')
    ch_idx = DATA.ALSUTRECHT.ica.icachansind;
else
    ch_idx = 1:size(DATA.data, 1);
end

[n_chans, n_pnts, n_trials] = size(DATA.data(ch_idx, :, :));
scalp_2d = double(reshape(DATA.data(ch_idx, :, :), n_chans, n_pnts * n_trials));

% Locate VEOG and HEOG channels
if ischar(cfg.veog_chan) || isstring(cfg.veog_chan)
    idx_veog = find(strcmpi({DATA.chanlocs.labels}, cfg.veog_chan), 1);
else
    idx_veog = cfg.veog_chan;
end

if ischar(cfg.heog_chan) || isstring(cfg.heog_chan)
    idx_heog = find(strcmpi({DATA.chanlocs.labels}, cfg.heog_chan), 1);
else
    idx_heog = cfg.heog_chan;
end

assert(~isempty(idx_veog), 'VEOG channel [%s] not found in DATA.', string(cfg.veog_chan));
assert(~isempty(idx_heog), 'HEOG channel [%s] not found in DATA.', string(cfg.heog_chan));

veog_2d = double(reshape(DATA.data(idx_veog, :, :), 1, n_pnts * n_trials));
heog_2d = double(reshape(DATA.data(idx_heog, :, :), 1, n_pnts * n_trials));

veog_3d = double(squeeze(DATA.data(idx_veog, :, :))); % [pnts x trials]
heog_3d = double(squeeze(DATA.data(idx_heog, :, :)));

% -------------------------------------------------------------------------
% 3. Global Continuous Correlation & Transmission Maps
% -------------------------------------------------------------------------
r_veog_cont    = zeros(n_chans, 1);
r_heog_cont    = zeros(n_chans, 1);
beta_veog_cont = zeros(n_chans, 1);
beta_heog_cont = zeros(n_chans, 1);

var_veog = var(veog_2d);
var_heog = var(heog_2d);

for c = 1:n_chans
    cov_v = cov(scalp_2d(c, :), veog_2d);
    r_mat_v = corrcoef(scalp_2d(c, :), veog_2d);
    r_veog_cont(c)    = r_mat_v(1, 2);
    beta_veog_cont(c) = cov_v(1, 2) / (var_veog + eps);

    cov_h = cov(scalp_2d(c, :), heog_2d);
    r_mat_h = corrcoef(scalp_2d(c, :), heog_2d);
    r_heog_cont(c)    = r_mat_h(1, 2);
    beta_heog_cont(c) = cov_h(1, 2) / (var_heog + eps);
end

% -------------------------------------------------------------------------
% 4. Dynamic IQR Event Gating (With Blink Precedence)
% -------------------------------------------------------------------------
% A. Vertical Blinks Detection
veog_detrend = veog_2d - movmedian(veog_2d, round(1.0 * fs));
veog_IQR     = iqr(veog_detrend);
veog_75P     = prctile(veog_detrend, 75);

thresh_veog   = max(cfg.min_phys_blink, veog_75P + (cfg.trIQRblink * veog_IQR));
blink_mask    = (veog_detrend > thresh_veog);
blink_mask_3d = reshape(blink_mask, n_pnts, n_trials);
blink_epochs  = any(blink_mask_3d, 1); % [1 x n_trials] boolean
has_blinks    = any(blink_epochs);

% B. Horizontal Saccades Detection (30 ms Haar Step Filter)
step_pts = max(2, round(0.030 * fs));
half_k   = ones(1, step_pts) / step_pts;
fwd_m    = conv(heog_2d, [half_k, zeros(1, step_pts)], 'same');
bwd_m    = conv(heog_2d, [zeros(1, step_pts), half_k], 'same');
heog_step = abs(fwd_m - bwd_m);

heog_IQR    = iqr(heog_step);
heog_75P    = prctile(heog_step, 75);

thresh_heog      = max(cfg.min_phys_saccade, heog_75P + (cfg.trIQRsaccade * heog_IQR));
saccade_mask_raw = (heog_step > thresh_heog);
saccade_mask_3d  = reshape(saccade_mask_raw, n_pnts, n_trials);

% Enforce Blink Precedence: Any epoch containing a blink is excluded from saccades
saccade_mask_3d(:, blink_epochs) = false;
saccade_mask   = saccade_mask_3d(:)';
saccade_epochs = any(saccade_mask_3d, 1);
has_saccades   = any(saccade_epochs);

% C. Event-Gated & Baseline Correlation Maps (Pre-Pruning)
r_veog_gated = zeros(n_chans, 1);
r_veog_base  = zeros(n_chans, 1);
delta_r_veog = zeros(n_chans, 1);

r_heog_gated = zeros(n_chans, 1);
r_heog_base  = zeros(n_chans, 1);
delta_r_heog = zeros(n_chans, 1);

if has_blinks && sum(blink_mask) > 10
    for c = 1:n_chans
        r_mat_e = corrcoef(scalp_2d(c, blink_mask), veog_2d(blink_mask));
        r_mat_b = corrcoef(scalp_2d(c, ~blink_mask), veog_2d(~blink_mask));
        r_veog_gated(c) = r_mat_e(1, 2);
        r_veog_base(c)  = r_mat_b(1, 2);
    end
    delta_r_veog = abs(r_veog_gated) - abs(r_veog_base);
else
    r_veog_gated = r_veog_cont;
    r_veog_base  = r_veog_cont;
end

if has_saccades && sum(saccade_mask) > 10
    for c = 1:n_chans
        r_mat_e = corrcoef(scalp_2d(c, saccade_mask), heog_2d(saccade_mask));
        r_mat_b = corrcoef(scalp_2d(c, ~saccade_mask), heog_2d(~saccade_mask));
        r_heog_gated(c) = r_mat_e(1, 2);
        r_heog_base(c)  = r_mat_b(1, 2);
    end
    delta_r_heog = abs(r_heog_gated) - abs(r_heog_base);
else
    r_heog_gated = r_heog_cont;
    r_heog_base  = r_heog_cont;
end

% -------------------------------------------------------------------------
% 5. Within-Epoch Event-Locked Snippet Averaging (ERP in uV)
% -------------------------------------------------------------------------
win_pts = round(cfg.win_sec * fs);
t_vec   = (-win_pts : win_pts) * (1000 / fs);

labels = {DATA.chanlocs(ch_idx).labels};

% Locate both left and right fronto-polar leads (BioSemi: C29 and C16)
idx_fp1 = find(strcmpi(labels, 'C29') | strcmpi(labels, 'Fp1') | strcmpi(labels, 'AF3'), 1);
idx_fp2 = find(strcmpi(labels, 'C16') | strcmpi(labels, 'Fp2') | strcmpi(labels, 'AF4'), 1);
if isempty(idx_fp1), idx_fp1 = 1; end
if isempty(idx_fp2), idx_fp2 = idx_fp1; end

idx_f7 = find(strcmpi(labels, 'D6') | strcmpi(labels, 'F7') | strcmpi(labels, 'AF7'), 1);
idx_f8 = find(strcmpi(labels, 'C6') | strcmpi(labels, 'F8') | strcmpi(labels, 'AF8'), 1);
if isempty(idx_f7), idx_f7 = 1; end
if isempty(idx_f8), idx_f8 = min(n_chans, 2); end

blink_snippets   = [];
saccade_snippets = [];

for tr = 1:n_trials
    v_tr = veog_3d(:, tr);
    h_tr = heog_3d(:, tr);

    % Blink Snippets
    if blink_epochs(tr)
        v_tr_det = v_tr - median(v_tr);
        [~, t_peak] = max(v_tr_det);
        if (t_peak - win_pts >= 1) && (t_peak + win_pts <= n_pnts)
            idx_w = (t_peak - win_pts) : (t_peak + win_pts);
            snip_v = double(DATA.data(ch_idx, idx_w, tr));
            snip_v = snip_v - mean(snip_v(:, [1:3, end-2:end]), 2);
            blink_snippets(:, :, end+1) = snip_v;
        end
    end

    % Saccade Snippets (Only evaluated on non-blink epochs)
    if saccade_epochs(tr) && ~blink_epochs(tr)
        fwd_tr = conv(h_tr, [half_k, zeros(1, step_pts)], 'same');
        bwd_tr = conv(h_tr, [zeros(1, step_pts), half_k], 'same');
        h_step_tr = abs(fwd_tr - bwd_tr);

        [~, t_step] = max(h_step_tr);
        if (t_step - win_pts >= 1) && (t_step + win_pts <= n_pnts)
            idx_w = (t_step - win_pts) : (t_step + win_pts);
            snip_h = double(DATA.data(ch_idx, idx_w, tr));
            snip_h = snip_h - mean(snip_h(:, [1:3, end-2:end]), 2);
            saccade_snippets(:, :, end+1) = snip_h;
        end
    end
end

if ~isempty(blink_snippets)
    erp_blink_scalp      = mean(blink_snippets, 3);
    max_blink_erp_peak   = max(abs(erp_blink_scalp(idx_fp1, :)));
    n_blink_snippets     = size(blink_snippets, 3);
else
    erp_blink_scalp      = zeros(n_chans, length(t_vec));
    max_blink_erp_peak   = 0;
    n_blink_snippets     = 0;
end

if ~isempty(saccade_snippets)
    erp_saccade_scalp    = mean(saccade_snippets, 3);
    lat_diff_erp         = erp_saccade_scalp(idx_f7, :) - erp_saccade_scalp(idx_f8, :);
    max_saccade_erp_peak = max(abs(lat_diff_erp));
    n_saccade_snippets   = size(saccade_snippets, 3);
else
    erp_saccade_scalp    = zeros(n_chans, length(t_vec));
    max_saccade_erp_peak = 0;
    n_saccade_snippets   = 0;
end

% -------------------------------------------------------------------------
% 6. Trial-by-Trial Residue Evaluation & Clean Epoch Definition
% -------------------------------------------------------------------------
trial_r_veog   = zeros(n_trials, 1);
trial_r_heog   = zeros(n_trials, 1);
trial_p2p_veog = zeros(n_trials, 1);
trial_p2p_heog = zeros(n_trials, 1);
blink_trials   = false(n_trials, 1);
saccade_trials = false(n_trials, 1);

for tr = 1:n_trials
    v_tr = veog_3d(:, tr);
    h_tr = heog_3d(:, tr);

    % fp1_tr   = double(DATA.data(ch_idx(idx_fp1), :, tr));
    lat_diff = double(DATA.data(ch_idx(idx_f7), :, tr) - DATA.data(ch_idx(idx_f8), :, tr));

    if blink_epochs(tr)
        blink_trials(tr) = true;

        % Check Left (Fp1) and Right (Fp2)
        r_fp1 = corrcoef(double(DATA.data(ch_idx(idx_fp1), :, tr)), v_tr);
        r_fp2 = corrcoef(double(DATA.data(ch_idx(idx_fp2), :, tr)), v_tr);

        p2p_fp1 = max(DATA.data(ch_idx(idx_fp1), :, tr)) - min(DATA.data(ch_idx(idx_fp1), :, tr));
        p2p_fp2 = max(DATA.data(ch_idx(idx_fp2), :, tr)) - min(DATA.data(ch_idx(idx_fp2), :, tr));

        % Store the worst-case metric between both eyes
        trial_r_veog(tr)   = max(r_fp1(1, 2), r_fp2(1, 2));
        trial_p2p_veog(tr) = max(p2p_fp1, p2p_fp2);
    end

    % Saccade check strictly on epochs without blinks
    if saccade_epochs(tr) && ~blink_epochs(tr)
        saccade_trials(tr) = true;
        r_tmp = corrcoef(lat_diff, h_tr);
        trial_r_heog(tr)   = r_tmp(1, 2);
        trial_p2p_heog(tr) = max(lat_diff) - min(lat_diff);
    end
end

% Adjusted trial thresholds
thresh_r_veog_trial = max(0.35, max(abs(r_veog_base)) + 0.15);
thresh_r_heog_trial = max(0.55, max(abs(r_heog_base)) + 0.15);

bad_blink_trials   = find(blink_trials   & (trial_r_veog > thresh_r_veog_trial) & (trial_p2p_veog > 15));
bad_saccade_trials = find(saccade_trials & (trial_r_heog > thresh_r_heog_trial) & (trial_p2p_heog > 15));

% All compromised trials
bad_trials_all = unique([bad_blink_trials(:); bad_saccade_trials(:)]);
clean_trials   = setdiff(1:n_trials, bad_trials_all);
n_clean_trials = length(clean_trials);

% -------------------------------------------------------------------------
% 7. Post-Pruning Clean Data Metrics
% -------------------------------------------------------------------------
r_veog_clean = zeros(n_chans, 1);
r_heog_clean = zeros(n_chans, 1);

if ~isempty(clean_trials)
    scalp_clean_2d = double(reshape(DATA.data(ch_idx, :, clean_trials), n_chans, n_pnts * n_clean_trials));
    veog_clean_2d  = double(reshape(DATA.data(idx_veog, :, clean_trials), 1, n_pnts * n_clean_trials));
    heog_clean_2d  = double(reshape(DATA.data(idx_heog, :, clean_trials), 1, n_pnts * n_clean_trials));

    for c = 1:n_chans
        r_v = corrcoef(scalp_clean_2d(c, :), veog_clean_2d);
        r_veog_clean(c) = r_v(1, 2);
        r_h = corrcoef(scalp_clean_2d(c, :), heog_clean_2d);
        r_heog_clean(c) = r_h(1, 2);
    end
end

% -------------------------------------------------------------------------
% 8. Pack Statistics Structure
% -------------------------------------------------------------------------
ocular_stats.r_veog_cont          = r_veog_cont;
ocular_stats.r_heog_cont          = r_heog_cont;
ocular_stats.r_veog_gated         = r_veog_gated;
ocular_stats.r_heog_gated         = r_heog_gated;
ocular_stats.r_veog_clean         = r_veog_clean;
ocular_stats.r_heog_clean         = r_heog_clean;
ocular_stats.delta_r_veog         = delta_r_veog;
ocular_stats.delta_r_heog         = delta_r_heog;
ocular_stats.beta_veog            = beta_veog_cont;
ocular_stats.beta_heog            = beta_heog_cont;
ocular_stats.max_abs_r_veog       = max(abs(r_veog_gated));
ocular_stats.max_abs_r_heog       = max(abs(r_heog_gated));
ocular_stats.max_abs_r_veog_clean = max(abs(r_veog_clean));
ocular_stats.max_abs_r_heog_clean = max(abs(r_heog_clean));
ocular_stats.max_delta_r_veog     = max(abs(delta_r_veog));
ocular_stats.max_delta_r_heog     = max(abs(delta_r_heog));
ocular_stats.max_blink_erp_peak   = max_blink_erp_peak;
ocular_stats.max_saccade_erp_peak = max_saccade_erp_peak;
ocular_stats.thresh_veog          = thresh_veog;
ocular_stats.thresh_heog          = thresh_heog;
ocular_stats.has_blinks           = has_blinks;
ocular_stats.has_saccades         = has_saccades;
ocular_stats.n_blink_snippets     = n_blink_snippets;
ocular_stats.n_saccade_snippets   = n_saccade_snippets;
ocular_stats.bad_blink_trials     = bad_blink_trials;
ocular_stats.bad_saccade_trials   = bad_saccade_trials;
ocular_stats.bad_trials_all       = bad_trials_all;
ocular_stats.clean_trials         = clean_trials;
ocular_stats.n_blink_trials       = sum(blink_trials);
ocular_stats.n_saccade_trials     = sum(saccade_trials);

% -------------------------------------------------------------------------
% 9. Console Diagnostic Summary
% -------------------------------------------------------------------------
fprintf('\n===================================================================\n');
fprintf('  Ocular Residual Diagnostic Assessment (Epoched Data)\n');
fprintf('===================================================================\n');
fprintf('  Total Epochs Analyzed: %d\n', n_trials);
fprintf('  Retained Clean Epochs: %d / %d (%.1f%%)\n', ...
    n_clean_trials, n_trials, (n_clean_trials / n_trials) * 100);
fprintf('  -----------------------------------------------------------------\n');
fprintf('  VERTICAL BLINKS (VEOG Threshold: %.1f uV):\n', thresh_veog);
fprintf('    - Active Epochs:                 %d (%.1f%% of recording time)\n', ...
    sum(blink_trials), (sum(blink_mask)/length(blink_mask))*100);
fprintf('    - Continuous Max |r|:            %.3f\n', max(abs(r_veog_cont)));
if has_blinks
    fprintf('    - Event-Gated Max |r|:           %.3f (Threshold: < %.2f)\n', ocular_stats.max_abs_r_veog, cfg.max_corr);
    fprintf('    - Contrast Delta-R (|r_e|-|r_b|): %.3f (Target: < 0.10)\n', ocular_stats.max_delta_r_veog);
    fprintf('    - Event-Locked ERP Peak (Fp1):   %.2f uV (Target: < 2.5 uV, N=%d)\n', max_blink_erp_peak, n_blink_snippets);
    fprintf('    - Compromised Blink Epochs:      %d / %d (%.1f%%)\n', ...
        length(bad_blink_trials), sum(blink_trials), (length(bad_blink_trials)/max(1, sum(blink_trials)))*100);
else
    fprintf('    - Event-Gated Assessment:        N/A (No macroscopic blinks detected)\n');
end

fprintf('  -----------------------------------------------------------------\n');
fprintf('  HORIZONTAL SACCADES (HEOG Step Threshold: %.1f uV | Non-Blink):\n', thresh_heog);
fprintf('    - Active Epochs:                 %d (%.1f%% of recording time)\n', ...
    sum(saccade_trials), (sum(saccade_mask)/length(saccade_mask))*100);
fprintf('    - Continuous Max |r|:            %.3f\n', max(abs(r_heog_cont)));
if has_saccades
    fprintf('    - Event-Gated Max |r|:           %.3f (Threshold: < %.2f)\n', ocular_stats.max_abs_r_heog, cfg.max_corr);
    fprintf('    - Contrast Delta-R (|r_e|-|r_b|): %.3f (Target: < 0.10)\n', ocular_stats.max_delta_r_heog);
    fprintf('    - Event-Locked ERP Peak (F7-F8): %.2f uV (Target: < 2.5 uV, N=%d)\n', max_saccade_erp_peak, n_saccade_snippets);
    fprintf('    - Compromised Saccade Epochs:    %d / %d (%.1f%%)\n', ...
        length(bad_saccade_trials), sum(saccade_trials), (length(bad_saccade_trials)/max(1, sum(saccade_trials)))*100);
else
    fprintf('    - Event-Gated Assessment:        N/A (No non-blink saccades detected)\n');
    fprintf('    - Diagnostic Note: Continuous R=%.2f reflects baseline cortical crosstalk into outer canthi leads.\n', max(abs(r_heog_cont)));
end
fprintf('  -----------------------------------------------------------------\n');
fprintf('  POST-PRUNING CLEAN EPOCHS (N = %d):\n', n_clean_trials);
fprintf('    - Clean VEOG Residual Max |r|:   %.3f\n', ocular_stats.max_abs_r_veog_clean);
fprintf('    - Clean HEOG Residual Max |r|:   %.3f\n', ocular_stats.max_abs_r_heog_clean);
fprintf('===================================================================\n\n');

% -------------------------------------------------------------------------
% 10. Diagnostic Figure Visualisation (2 Rows x 2 Columns)
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot, return; end

fh = figure('Color', 'w', 'Position', [100, 100, 1050, 800], 'Name', 'Ocular Residual Diagnostics (Pre vs Post Pruning)');

% Row 1, Col 1: Pre-Pruning VEOG Topography (Delta-R or Gated)
sbh1 = subplot(2, 2, 1);
if has_blinks
    title_v1 = sprintf('Pre-Pruning VEOG Delta-R (ERP = %.1f uV)', max_blink_erp_peak);
    render_topoplot(delta_r_veog, DATA.chanlocs(ch_idx), title_v1, sbh1, [-0.20, 0.20]);
else
    title_v1 = sprintf('Pre-Pruning VEOG (|r|_{max} = %.2f, No Blinks)', max(abs(r_veog_cont)));
    render_topoplot(r_veog_cont, DATA.chanlocs(ch_idx), title_v1, sbh1, [-0.20, 0.20]);
end

% Row 1, Col 2: Pre-Pruning HEOG Topography (Delta-R or Gated)
sbh2 = subplot(2, 2, 2);
if has_saccades
    title_h1 = sprintf('Pre-Pruning HEOG Delta-R (ERP = %.1f uV)', max_saccade_erp_peak);
    render_topoplot(delta_r_heog, DATA.chanlocs(ch_idx), title_h1, sbh2, [-0.20, 0.20]);
else
    title_h1 = sprintf('Pre-Pruning HEOG Crosstalk (|r|_{max} = %.2f)', max(abs(r_heog_cont)));
    render_topoplot(r_heog_cont, DATA.chanlocs(ch_idx), title_h1, sbh2, [-0.20, 0.20]);
end

% Row 2, Col 1: Post-Pruning Clean VEOG Topography
sbh3 = subplot(2, 2, 3);
title_v2 = sprintf('Post-Pruning Clean VEOG (N = %d Epochs, |r|_{max} = %.2f)', ...
    n_clean_trials, ocular_stats.max_abs_r_veog_clean);
render_topoplot(r_veog_clean, DATA.chanlocs(ch_idx), title_v2, sbh3, [-0.20, 0.20]);

% Row 2, Col 2: Post-Pruning Clean HEOG Topography
sbh4 = subplot(2, 2, 4);
title_h2 = sprintf('Post-Pruning Clean HEOG (N = %d Epochs, |r|_{max} = %.2f)', ...
    n_clean_trials, ocular_stats.max_abs_r_heog_clean);
render_topoplot(r_heog_clean, DATA.chanlocs(ch_idx), title_h2, sbh4, [-0.20, 0.20]);

drawnow;

end

% =========================================================================
% Local Helper: Safe Topoplot Dispatcher
% =========================================================================
function render_topoplot(values, chanlocs, plot_title, ax_handle, map_limits)
axes(ax_handle);
if exist('mytopoplot', 'file') == 2
    mytopoplot(values, [], plot_title, ax_handle, map_limits);
    colorbar;
elseif exist('topoplot', 'file') == 2
    topoplot(values, chanlocs, 'maplimits', map_limits, 'electrodes', 'off', 'style', 'map');
    title(plot_title, 'FontSize', 10, 'FontWeight', 'bold');
    colorbar;
else
    bar(values);
    title(plot_title, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Correlation Value');
    xlabel('Channels');
end
end