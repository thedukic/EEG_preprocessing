function EEG = do_dss_blink(EEG, cfg)
% DO_DSS_BLINK Detects ocular excursions on continuous 2D EEG using the dedicated
% VEOG channel, extracts blink artefact subspace via Denoising Source Separation (DSS),
% regresses blink components out of EEG channels, and renders topographic verification.

if nargin < 2, cfg = struct(); end

% 1. Validation & Default Configuration
assert(ismatrix(EEG.data), 'EEG.data must be 2D continuous data. Run this function prior to epoching.');

if ~isfield(cfg, 'veog_chan'),   cfg.veog_chan   = 'VEOG'; end
if ~isfield(cfg, 'iqr_mult'),    cfg.iqr_mult    = 3.0;    end
if ~isfield(cfg, 'min_deflect'), cfg.min_deflect = 50;     end
if ~isfield(cfg, 'pad_sec'),     cfg.pad_sec     = 0.25;   end
if ~isfield(cfg, 'n_remove'),    cfg.n_remove    = [];     end % Explicit override (e.g., 1 or 2)
if ~isfield(cfg, 'do_plot'),     cfg.do_plot     = true;   end

fs = EEG.srate;

% 2. VEOG Channel Extraction
if ischar(cfg.veog_chan) || isstring(cfg.veog_chan)
    idx_veog = find(strcmpi({EEG.chanlocs.labels}, cfg.veog_chan), 1);
else
    idx_veog = cfg.veog_chan;
end
assert(~isempty(idx_veog), 'VEOG channel [%s] not found in EEG.', string(cfg.veog_chan));

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

fprintf('--- Continuous DSS Blink Cleaning ---\n');
fprintf('  Detector Channel: %s\n', cfg.veog_chan);
fprintf('  Threshold:        %.2f uV\n', det_thresh);
fprintf('  Artefact masked:  %d samples (%.2f s)\n', n_art_samples, n_art_samples / fs);
fprintf('  Clean retained:   %d samples (%.2f s)\n', n_clean_samples, n_clean_samples / fs);

% Initialise logging struct
EEG.ALSUTRECHT.dss.status      = false;
EEG.ALSUTRECHT.dss.n_removed   = 0;
EEG.ALSUTRECHT.dss.removed_idx = [];

% Minimum data checks
min_req_samples = round(2.0 * fs);
if n_art_samples < min_req_samples
    warning('Insufficient blink data detected (< 2.0 s). Skipping DSS.');
    return;
end
if n_clean_samples < min_req_samples
    warning('Insufficient clean baseline segments remaining. Skipping DSS.');
    return;
end

% 4. DSS Covariance Estimation & Component Decomposition
mask_eeg = strcmpi({EEG.chanlocs.type}, 'EEG');
raw_eeg  = double(EEG.data(mask_eeg, :));

% Sample 30s to determine numerical rank (guards spherical whitening in nt_dss0)
data_rank = rank(raw_eeg);
data_rank = min(data_rank, 20);

% 1. Create a temporary low-pass filtered copy strictly for covariance estimation
% (e.g., 4th-order zero-phase Butterworth at 5 Hz)
[b_lp, a_lp] = butter(4, 5 / (fs / 2), 'low');
eeg_low      = filtfilt(b_lp, a_lp, raw_eeg')';

% 2. Estimate C1 and C0 using the low-passed data
C1 = nt_cov(bsxfun(@times, eeg_low', double(mask_blink')));
C0 = nt_cov(bsxfun(@times, eeg_low', double(~mask_blink')));

% Compute DSS unmixing matrix keeping up to active data rank
[todss, pwr0, pwr1] = nt_dss0(C0, C1, data_rank, []);
p1 = (pwr1 ./ pwr0) / sum(pwr1 ./ pwr0);
fromdss = pinv(todss);

% Transform scalp data to DSS component activations
noise_components = todss' * raw_eeg; % [n_comp x n_samples]

% 1. Define anterior and posterior channel zones from montage coordinates
templates_ica = load_ictemplateweights(EEG);
artifact_templates = [templates_ica.Blinkweights0, templates_ica.Blinkweights1];

num_check = min(12, size(fromdss, 1));

cfg_tmp.max_check  = num_check;
cfg_tmp.blink      = cfg.roi.blink;
cfg_tmp.blink_anti = cfg.roi.blink_anti;
[indx_remove, sim_scores, ap_scores, report] = match_dss_blink(fromdss, EEG.chanlocs, artifact_templates, cfg_tmp);

% Ensure components are strictly within the top 5 (blinks are low-rank)
indx_remove = indx_remove(indx_remove <= 3);

% -------------------------------------------------------------------------
% Render Topoplots
% -------------------------------------------------------------------------
if isempty(indx_remove)
    warning('No DSS blink components identified. Skipping cleaning.');
    fh = figure('Color', 'w', 'Position', [100, 100, 1400, 340], 'Name', 'DSS Blink Power Verification');
    tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i_cmp = 1:5
        if ismember(i_cmp, indx_remove)
            tag = '*';
        else
            tag = '';
        end

        title_topo = sprintf('%sDSS%d\nR = %.2f, AP = %.2f%s', ...
            tag, i_cmp, sim_scores(i_cmp), ap_scores(i_cmp), tag);

        mytopoplot(fromdss(i_cmp, :), [], title_topo, nexttile);
    end

    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_dss_blink_power'], [36 10]);
    return;
end

fprintf('  DSS components flagged for removal: %s\n', mat2str(indx_remove(:)'));

% 5. Artefact Subspace Regression
eeg_clean = nt_tsr(raw_eeg', noise_components(indx_remove, :)')';

% 6. Topographic Power Diagnostics
if cfg.do_plot
    raw_art   = raw_eeg(:, mask_blink);
    clean_art = eeg_clean(:, mask_blink);
    raw_art   = raw_art - mean(raw_art, 2);
    clean_art = clean_art - mean(clean_art, 2);

    pow_pre  = mean(raw_art.^2, 2);
    pow_post = mean(clean_art.^2, 2);
    delta_db_blink = 10 * log10(pow_post ./ pow_pre);

    raw_clean_seg   = raw_eeg(:, ~mask_blink);
    clean_clean_seg = eeg_clean(:, ~mask_blink);
    raw_clean_seg   = raw_clean_seg - mean(raw_clean_seg, 2);
    clean_clean_seg = clean_clean_seg - mean(clean_clean_seg, 2);

    pow_clean_pre  = mean(raw_clean_seg.^2, 2);
    pow_clean_post = mean(clean_clean_seg.^2, 2);
    delta_db_clean = 10 * log10(pow_clean_post ./ pow_clean_pre);

    fh = figure('Color', 'w', 'Position', [100, 100, 1400, 340], 'Name', 'DSS Blink Power Verification');
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    max_raw_pow = max(pow_pre);

    % 1. Pre-Cleaning Blink Power
    ax1 = nexttile;
    mytopoplot(pow_pre, [], '', ax1, [0, max_raw_pow]);
    hcb1 = colorbar(ax1);
    hcb1.Title.String = '\muV^2';
    title(ax1, sprintf('Pre-DSS Blink Power\n(Max: %.1f \\muV^2)', max_raw_pow), 'FontWeight', 'bold');

    % 2. Post-Cleaning Blink Power
    ax2 = nexttile;
    mytopoplot(pow_post, [], '', ax2, [0, max_raw_pow]);
    hcb2 = colorbar(ax2);
    hcb2.Title.String = '\muV^2';
    title(ax2, sprintf('Post-DSS Blink Power\n(Max: %.1f \\muV^2)', max(pow_post)), 'FontWeight', 'bold');

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

    for i_cmp = 1:4
        if ismember(i_cmp, indx_remove)
            col = '\color{red}';
        else
            col = '\color{black}';
        end

        title_topo = sprintf('%sDSS%d\n%sR = %.2f, AP = %.2f', ...
            col, i_cmp, col, sim_scores(i_cmp), ap_scores(i_cmp));

        mytopoplot(fromdss(i_cmp, :), [], title_topo, nexttile);
    end

    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_dss_blink_power'], [36 12]);
end

% 7. Validate Output & Store
flag_successful = true;
if ~isreal(eeg_clean) || any(isnan(eeg_clean(:)))
    warning('Something went wrong with DSS: output contains NaNs or complex values.');
    flag_successful = false;
end

if flag_successful
    EEG.data(mask_eeg, :) = eeg_clean;

    EEG.ALSUTRECHT.dss.status      = true;
    EEG.ALSUTRECHT.dss.n_removed   = length(indx_remove);
    EEG.ALSUTRECHT.dss.removed_idx = indx_remove;
    EEG.ALSUTRECHT.dss.scores      = p1;
    EEG.ALSUTRECHT.dss.report      = report;
    fprintf('DSS Blink Cleaning Completed: %d component(s) removed.\n', length(indx_remove));
end

end


function [bad_cmps, sim_scores, ap_scores, report] = match_dss_blink(fromdss, chanlocs, templates, cfg)
% MATCH_DSS_BLINK Identifies blink components from DSS spatial filters (BioSemi 128)
%
% Syntax:
%   [bad_cmps, sim_scores, ap_scores, report] = match_dss_blink(fromdss, chanlocs, templates, cfg)
%
% Inputs:
%   fromdss   : [num_cmps x 128] or [128 x num_cmps] DSS mixing matrix (pinv(todss))
%   chanlocs  : EEGLAB chanlocs struct array (EEG.chanlocs)
%   templates : Struct containing blink weights (e.g. templatesICA.Blinkweights0)
%               OR numeric matrix of template columns [128 x num_variants]
%   cfg       : (Optional) Configuration structure:
%               - cfg.sim_thresh : Minimum spatial similarity (default: 0.65)
%               - cfg.ap_thresh  : Minimum antero-posterior ratio (default: 1.25)
%               - cfg.max_check  : Maximum component index to evaluate (default: 5)
%               - cfg.metric     : 'cosine' (default) or 'pearson'
%
% Outputs:
%   bad_cmps   : Row vector of identified blink component indices (e.g. [1, 2])
%   sim_scores : [1 x num_checked] Best similarity score per component
%   ap_scores  : [1 x num_checked] Antero-posterior voltage ratio per component
%   report     : Summary structure with detailed metrics per component

if nargin < 4, cfg = struct(); end

% 1. Default Parameters
if ~isfield(cfg, 'sim_thresh'), cfg.sim_thresh = 0.80;     end
if ~isfield(cfg, 'ap_thresh'),  cfg.ap_thresh  = 1.25;     end
if ~isfield(cfg, 'max_check'),  cfg.max_check  = 10;        end
if ~isfield(cfg, 'metric'),     cfg.metric     = 'cosine'; end

% 2. Orient Spatial Patterns to [128 x num_cmps]
fromdss = double(fromdss);
if size(fromdss, 1) < size(fromdss, 2)
    % Transpose from [num_cmps x 128] to [128 x num_cmps]
    W_cmp = fromdss';
else
    W_cmp = fromdss;
end

num_chans = size(W_cmp, 1);
num_cmps  = size(W_cmp, 2);
num_check = min(cfg.max_check, num_cmps);
W_eval    = W_cmp(:, 1:num_check);

% 3. Extract Blink Template Variants
if isstruct(templates)
    fnames = fieldnames(templates);
    blink_fields = fnames(contains(lower(fnames), 'blink'));
    assert(~isempty(blink_fields), 'No fields containing "blink" found in templates struct.');

    W_tgt = [];
    for i_f = 1:length(blink_fields)
        tmpl = double(templates.(blink_fields{i_f}));
        assert(size(tmpl, 1) == num_chans, ...
            'Template %s has %d channels; expected %d.', blink_fields{i_f}, size(tmpl, 1), num_chans);
        W_tgt = [W_tgt, tmpl]; %#ok<AGROW>
    end
else
    W_tgt = double(templates);
    assert(size(W_tgt, 1) == num_chans, ...
        'Input template matrix has %d channels; expected %d.', size(W_tgt, 1), num_chans);
end

% 4. Multi-Variant Spatial Similarity
switch lower(cfg.metric)
    case 'cosine'
        % Normalised dot product: |u' * v| / (||u|| * ||v||)
        norm_cmp = sqrt(sum(W_eval.^2, 1));
        norm_tgt = sqrt(sum(W_tgt.^2, 1))';
        dot_prod = abs(W_tgt' * W_eval);
        sim_mat  = dot_prod ./ (norm_tgt * norm_cmp + eps);
    case 'pearson'
        sim_mat  = abs(corr(W_eval, W_tgt, 'Type', 'Pearson'))';
    otherwise
        error('Unrecognised metric "%s". Use ''cosine'' or ''pearson''.', cfg.metric);
end
[sim_scores, best_variant] = max(sim_mat, [], 1);

% 5. Anterior Dominance & Occipital Rejection Gate
chan_labels    = upper({chanlocs(1:num_chans).labels});
mask_anterior  = ismember(chan_labels, cfg.blink);
mask_posterior = ismember(chan_labels, cfg.blink_anti);

ap_scores        = zeros(1, num_check);
is_peak_anterior = false(1, num_check);

for i_c = 1:num_check
    topo = W_eval(:, i_c);

    % Peak absolute deflection must lie in anterior zone
    % Do this on ALL channels that are anteriror, X > 0
    [~, peak_idx] = max(abs(topo));
    is_peak_anterior(i_c) = mask_anterior(peak_idx);

    % Antero-posterior voltage amplitude ratio
    ap_scores(i_c) = max(abs(topo(mask_anterior))) / (max(abs(topo(mask_posterior))));
end

% 6. Joint Classification
pass_sim      = sim_scores >= cfg.sim_thresh;
% pass_spatial  = is_peak_anterior & (ap_scores >= cfg.ap_thresh);
pass_spatial  = ap_scores >= cfg.ap_thresh;
bad_mask      = pass_sim & pass_spatial;
bad_cmps      = find(bad_mask);

% 7. Diagnostic Packaging
report                 = struct();
report.bad_cmps        = bad_cmps;
report.sim_scores      = sim_scores;
report.ap_scores       = ap_scores;
report.best_variant    = best_variant;
report.is_peak_frontal = is_peak_anterior;

if isempty(bad_cmps)
    fprintf('[DSS Blink Match] Screened top %d components: identified none (Sim >= %.2f, AP >= %.2f)\n', num_check, cfg.sim_thresh, cfg.ap_thresh);
    fprintf('[DSS Blink Match: Fallback] Testing Component 1 for atypical/unilateral blink...\n');

    w_cand = W_eval(:, 1);
    active_labels = {chanlocs.labels};

    % Extract predefined blink channels dynamically and split by Y coordinate (+Y = Left, -Y = Right)
    blink_mask = ismember(active_labels, cfg.blink);
    blink_locs = chanlocs(blink_mask);

    fp_left_names  = {blink_locs([blink_locs.Y] > 0).labels};
    fp_right_names = {blink_locs([blink_locs.Y] < 0).labels};

    idx_l = find(ismember(active_labels, fp_left_names));
    idx_r = find(ismember(active_labels, fp_right_names));

    % Posterior channels defined dynamically by coordinates (posterior to vertex: X < 0)
    % idx_post = [chanlocs.X] < 0;
    idx_post = ismember(active_labels, cfg.blink_anti);

    % Mean RMS per eye cluster
    rms_left  = sqrt(mean(w_cand(idx_l).^2));
    rms_right = sqrt(mean(w_cand(idx_r).^2));
    rms_post  = sqrt(mean(w_cand(idx_post).^2)) + eps;

    rms_eye_max = max(rms_left, rms_right);
    eye_to_post_ratio = rms_eye_max / rms_post;

    % Check unipolar sign across whichever eye cluster is dominant
    if rms_left >= rms_right
        active_vals = w_cand(idx_l);
    else
        active_vals = w_cand(idx_r);
    end
    is_same_sign = all(active_vals > 0) || all(active_vals < 0);

    % Total fronto-polar energy fraction from cfg.blink
    total_energy = sum(w_cand.^2) + eps;
    fp_energy    = sum(w_cand([idx_l, idx_r]).^2);
    fp_fraction  = fp_energy / total_energy;

    % Decision gate
    if (eye_to_post_ratio > 3.5) && is_same_sign && (fp_fraction > 0.35)
        fprintf('  Confirmed asymmetric blink in Component 1 (Eye/Post ratio: %.1fx, FP energy: %.1f%%)\n', ...
            eye_to_post_ratio, fp_fraction * 100);
        bad_cmps = 1;
    else
        fprintf('  Component 1 rejected (Eye/Post ratio: %.1fx, Same sign: %d, FP energy: %.1f%%)\n', ...
            eye_to_post_ratio, is_same_sign, fp_fraction * 100);
    end

    % Log
    report.dss1.eye_to_post_ratio    = eye_to_post_ratio;
    report.dss1.is_same_sign      = is_same_sign;
    report.dss1.fp_fraction       = fp_fraction;
else
    fprintf('[DSS Blink Match] Screened top %d components: identified %s (Sim >= %.2f, AP >= %.2f)\n', ...
        num_check, mat2str(bad_cmps), cfg.sim_thresh, cfg.ap_thresh);
end

end