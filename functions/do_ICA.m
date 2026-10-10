function EEG = do_ica(EEG, cfg)
% DO_ICA
% Runs ICA (RunICA, Picard, FastICA, or AMICA) on scalp EEG channels using a
% temporary 1.0 Hz high-pass copy for decomposition training, transferring the
% spatial weights back to the base dataset (preserving 0.1 Hz delta & 1/f slope).

fprintf('\n==================================================\n');
fprintf('ICA: Independent Component Analysis (Two-Track 1.0 Hz)\n');
fprintf('==================================================\n');

% -------------------------------------------------------------------------
% 1. Channel selection
% -------------------------------------------------------------------------
eeg_chan_idx = find(strcmpi({EEG.chanlocs.type}, 'EEG'));

% -------------------------------------------------------------------------
% 2. Create Temporary 1.0 Hz Dataset for ICA Training
% -------------------------------------------------------------------------
fprintf('Preparing temporary 1.0 Hz filtered dataset for ICA training...\n');

cfg_tmp = cfg;
cfg_tmp.flt.erp.hp = [1 4];
cfg_tmp.flt.rs.hp  = [1 4];
cfg_tmp.flt.mt.hp  = [1 4];

EEG_train = do_filtering(EEG, 'highpass', 'eeg', cfg_tmp);

eeg_data_2d = double(reshape(EEG_train.data(eeg_chan_idx, :, :), length(eeg_chan_idx), []));
[n_chan, n_pnts] = size(eeg_data_2d);

% -------------------------------------------------------------------------
% 3. Spatial Covariance, SVD & PCA Rank Verification (on 1.0 Hz data)
% -------------------------------------------------------------------------
cov_mat = (eeg_data_2d * eeg_data_2d') / (n_pnts - 1);
[V, D]  = eig(cov_mat);
[eig_vals, sort_idx] = sort(diag(D), 'descend');
V = V(:, sort_idx);

explained    = eig_vals ./ sum(eig_vals);
explainedTmp = cumsum(explained);
num_pca      = find(explainedTmp >= 0.95, 1);

% Determine data rank
num_rank = get_rank(eeg_data_2d);
[num_ica, num_ica_max, k_achieved] = select_num_ics(EEG_train, cfg.ica.num_ica);

if num_ica > num_rank
    fprintf('Warning: Requested %d ICs exceeds data rank of %d. Adjusting to rank.\n', num_ica, num_rank);
    num_ica = num_rank;
end

var_pca = 100 * sum(explained(1:num_ica));

% -------------------------------------------------------------------------
% 4. Plot Kept vs Discarded PCA diagnostics
% -------------------------------------------------------------------------
if num_ica < n_chan
    fh = plot_pca_residuals(EEG_train, eeg_chan_idx, eeg_data_2d, V, eig_vals, num_ica, 8);
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_power_removed'], [20 20]);
end

% -------------------------------------------------------------------------
% 5. Execute Decomposition on 1.0 Hz Data
% -------------------------------------------------------------------------
if strcmpi(cfg.ica.type, 'AMICA')
    if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'subject') && isfield(EEG.ALSUTRECHT.subject, 'icadata')
        outdir = EEG.ALSUTRECHT.subject.icadata;
    else
        outdir = fullfile(pwd, 'amica_output');
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    fprintf('Running AMICA (outdir: %s)...\n', outdir);

    runamica15(eeg_data_2d, 'num_chans', n_chan, 'outdir', outdir, ...
        'pcakeep', num_ica, 'num_models', 1, 'numprocs', 1, 'max_threads', 10, ...
        'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1, 'max_iter', 1000);

    modout        = loadmodout15(outdir);
    icaweights    = modout.W;
    icasphere     = modout.S(1:modout.num_pcs, :);
    icawinv       = modout.A;
    EEG.etc.amica = modout;
else
    % Standard pop_runica on the 1.0 Hz training set
    EEG_train = pop_runica(EEG_train, 'icatype', lower(cfg.ica.type), 'extended', 1, ...
        'pca', num_ica, 'chanind', eeg_chan_idx, 'lrate', 1e-4, 'maxsteps', 2500);

    icaweights = EEG_train.icaweights;
    icasphere  = EEG_train.icasphere;
    icawinv    = EEG_train.icawinv;
end

% Clear the temporary 1.0 Hz dataset from memory
clearvars EEG_train;

% -------------------------------------------------------------------------
% 6. Transfer Decomposition Matrices Back to Base EEG (0.1 Hz)
% -------------------------------------------------------------------------
EEG.icaweights  = icaweights;
EEG.icasphere   = icasphere;
EEG.icawinv     = icawinv;
EEG.icachansind = eeg_chan_idx;

% Validate structure and clear cached activations to conserve memory
EEG = eeg_checkset(EEG, 'ica');
EEG.icaact = [];

% -------------------------------------------------------------------------
% 7. Post-ICA Variance Tracking (Computed on the Preserved Base Data)
% -------------------------------------------------------------------------
[vaf_compvar, var_global, var_peak_chan, peak_chan_idx] = calc_ic_variance(EEG);

% Cache decomposition details
EEG.ALSUTRECHT.ica.icasphere       = EEG.icasphere;
EEG.ALSUTRECHT.ica.icaweights      = EEG.icaweights;
EEG.ALSUTRECHT.ica.icawinv         = EEG.icawinv;
EEG.ALSUTRECHT.ica.icachansind     = EEG.icachansind;
EEG.ALSUTRECHT.ica.num_req         = cfg.ica.num_ica;
EEG.ALSUTRECHT.ica.num_max         = num_ica_max;
EEG.ALSUTRECHT.ica.num_done        = num_ica;
EEG.ALSUTRECHT.ica.num_pca         = num_pca;
EEG.ALSUTRECHT.ica.k_achieved      = k_achieved;
EEG.ALSUTRECHT.ica.var_pca         = var_pca;
EEG.ALSUTRECHT.ica.vaf_compvar     = vaf_compvar;
EEG.ALSUTRECHT.ica.var_ica         = var_global;
EEG.ALSUTRECHT.ica.var_peak_chan   = var_peak_chan;
EEG.ALSUTRECHT.ica.peak_chan_idx   = peak_chan_idx;

% Console and File Logging
log_text = sprintf(['\nICA Summary:\n', ...
    '  Algorithm:             %s\n', ...
    '  Training High-Pass:    1.0 Hz (Weights back-projected to base data)\n', ...
    '  Channels Analysed:     %d\n', ...
    '  Requested Max ICs:     %d\n', ...
    '  Ideal Max ICs:         %d\n', ...
    '  Estimated ICs:         %d\n', ...
    '  Achieved k:            %0.1f\n', ...
    '  Variance Retained:     %.2f%%\n', ...
    '  95%% Variance Rank:     %d\n'], ...
    cfg.ica.type, length(eeg_chan_idx), max(cfg.ica.num_ica), num_ica_max, num_ica, k_achieved, var_pca, num_pca);

fprintf('%s', log_text);

if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'subject') && ...
        isfield(EEG.ALSUTRECHT.subject, 'fid') && ~isempty(EEG.ALSUTRECHT.subject.fid) && ...
        EEG.ALSUTRECHT.subject.fid > 0
    fprintf(EEG.ALSUTRECHT.subject.fid, '%s', log_text);
end

fprintf('\nDone!\n');

end

% =========================================================================
% Helper: Diagnostic Plot for Discarded PCA Subspace
% =========================================================================
function fh = plot_pca_residuals(EEG, eeg_chan_idx, eeg_data_2d, V, eig_vals, num_ica, n_discard_plot)

if nargin < 7 || isempty(n_discard_plot)
    n_discard_plot = 5;
end

n_chan = length(eeg_chan_idx);
fs     = EEG.srate;

% 1. Reconstruct Sensor Data from Kept and Discarded Subspaces
V_kept = V(:, 1:num_ica);
V_disc = V(:, (num_ica + 1):end);

data_kept = V_kept * (V_kept' * eeg_data_2d);
data_disc = V_disc * (V_disc' * eeg_data_2d);

% 2. Calculate Welch PSD
nfft     = min(size(eeg_data_2d, 2), 256 * 4);
window   = hanning(min(size(eeg_data_2d, 2), 256 * 2));
noverlap = floor(length(window) / 2);

[psd_orig_all, freqs] = pwelch(eeg_data_2d', window, noverlap, nfft, fs);
psd_kept_all          = pwelch(data_kept',   window, noverlap, nfft, fs);
psd_disc_all          = pwelch(data_disc',   window, noverlap, nfft, fs);

global_psd_orig = mean(psd_orig_all, 2);
global_psd_kept = mean(psd_kept_all, 2);
global_psd_disc = mean(psd_disc_all, 2);

f_max    = 70;
f_mask   = freqs >= 1  & freqs <= f_max;
f_broad  = freqs >= 2  & freqs <= 48;
beta_idx = freqs >= 13 & freqs <= 35;

pct_disc_power = (sum(eig_vals((num_ica + 1):end)) / sum(eig_vals)) * 100;

% 3. Relative Beta Power
rel_beta_kept = (sum(psd_kept_all(beta_idx, :), 1) ./ sum(psd_kept_all(f_broad, :), 1)) * 100;
rel_beta_disc = (sum(psd_disc_all(beta_idx, :), 1) ./ sum(psd_disc_all(f_broad, :), 1)) * 100;

% 4. Generate Diagnostic Figure
fh = figure('Color', 'w');

subplot(3, 4, [1, 2]);
plot(freqs(f_mask), 10*log10(global_psd_orig(f_mask)), 'k', 'LineWidth', 1.8); hold on;
plot(freqs(f_mask), 10*log10(global_psd_kept(f_mask)), 'b--', 'LineWidth', 1.4);
plot(freqs(f_mask), 10*log10(global_psd_disc(f_mask)), 'r', 'LineWidth', 1.5);
grid on;
xlim([1, f_max]);
xlabel('Frequency (Hz)');
ylabel('PSD (10*log10 \muV^2/Hz)');
title({'Global Power Spectrum', ...
    sprintf('Discarded ICs: %d-%d (%.2f%% Total power)', num_ica + 1, n_chan, pct_disc_power)});
legend({'Original data', 'Kept PCs', 'Discarded PCs'}, 'Location', 'northeast');
pbaspect([1.6 1 1]);

sph = subplot(3, 4, 3);
max_kept = max(rel_beta_kept);
min_kept = min(rel_beta_kept);
if max_kept - min_kept < 5
    min_kept = max_kept - 5;
end
mytopoplot(rel_beta_kept(:), [], {'Kept', 'Rel. \beta (13-35 Hz)'}, sph, [min_kept, max_kept]);
cb = colorbar;
cb.Label.String   = '%';
cb.Location       = 'southoutside';
cb.Label.FontSize = 8;

scale_factor = 1.3;
pos = get(sph, 'Position');
new_w = pos(3) * scale_factor;
new_h = pos(4) * scale_factor;
new_left = pos(1) - (new_w - pos(3));
new_bottom = pos(2) - (new_h - pos(4));
set(sph, 'Position', [new_left, new_bottom, new_w, new_h]);

sph = subplot(3, 4, 4);
max_disc = max(rel_beta_disc);
min_disc = min(rel_beta_disc);
if max_disc - min_disc < 5
    min_disc = max_disc - 5;
end
mytopoplot(rel_beta_disc(:), [], {'Discarded', 'Rel. \beta (13-35 Hz)'}, sph, [min_disc, max_disc]);
cb = colorbar;
cb.Label.String   = '%';
cb.Location       = 'southoutside';
cb.Label.FontSize = 8;

pos = get(sph, 'Position');
new_w = pos(3) * scale_factor;
new_h = pos(4) * scale_factor;
new_left = pos(1) - (new_w - pos(3));
new_bottom = pos(2) - (new_h - pos(4));
set(sph, 'Position', [new_left, new_bottom, new_w, new_h]);

pcs_to_plot = (num_ica + 1) : min(num_ica + n_discard_plot, n_chan);

for i = 1:length(pcs_to_plot)
    pc_idx = pcs_to_plot(i);
    sph = subplot(3, 4, 4 + i);

    pc_weights = V(:, pc_idx);
    cmax = max(abs(pc_weights));
    if cmax == 0, cmax = 1; end

    var_exp   = (eig_vals(pc_idx) / sum(eig_vals)) * 100;
    title_str = sprintf('PC %d (Var: %.3f%%)', pc_idx, var_exp);

    mytopoplot(pc_weights, [], title_str, sph, [-cmax, cmax]);
end
end

% =========================================================================
% Helper: Dynamic Selection of Components
% =========================================================================
function [num_ica, num_ica_max, k_achieved] = select_num_ics(EEG, num_ica_candidates)
% SELECT_NUM_ICS
% Dynamically selects the number of PCs for ICA based on the SCCN/Makoto heuristic:
%   Data Points Needed = 20 * (N^2)  -->  N <= sqrt(Data Points / 20)
%
% Handles scalar inputs (e.g., 80) or candidate vectors (e.g., [128 100 80 70 50]).
% Caps selection by the recommended maximum, data rank, and available scalp channels.
%
% Paradigm support at 256 Hz (k = 20):
%   MMN   3*7     ~ 21 min      -> Supports up to 127 ICs (Select: 100-128)
%   SART  3*5     ~ 15 min      -> Supports up to 107 ICs (Select: 100)
%   RS    2x3x2   ~ 12 min      -> Supports up to 96 ICs  (Select: 80-100)
%   MT    7+3(+7) ~ 10 (17) min -> Supports up to 88 (114) ICs (Select: 80-100)
%   DUB   EO      ~ 6 min       -> Supports up to 68 ICs  (Select: 64-70)
%
% Minimum EEG needed (at 256 Hz, k = 20):
%   128 elecs: 128^2 * 20 / 256 / 60 ~ 21.3 min (327,680 pts)
%   100 PCs:   100^2 * 20 / 256 / 60 ~ 13.0 min (200,000 pts)
%   80 PCs:     80^2 * 20 / 256 / 60 ~  8.3 min (128,000 pts)
%   70 PCs:     70^2 * 20 / 256 / 60 ~  6.4 min  (98,000 pts)
%   50 PCs:     50^2 * 20 / 256 / 60 ~  3.3 min  (50,000 pts)
%   32 PCs:     32^2 * 20 / 256 / 60 ~  1.3 min  (20,480 pts)
%
% Reference: Makeig & Onton (2011), Oxford Handbook of ERP Components.
% See also:  https://eeglab.ucsd.edu/wiki/Makoto's_preprocessing_pipeline

if nargin < 2 || isempty(num_ica_candidates)
    num_ica_candidates = [128, 100, 80, 70, 64, 50, 32];
end

fprintf('Selecting #PCs for ICA (SCCN Heuristic k = 20)...\n');
fprintf('Requested max candidate: %d\n', max(num_ica_candidates));

if isfield(EEG, 'pnts') && isfield(EEG, 'trials') && ~isempty(EEG.pnts) && ~isempty(EEG.trials)
    data_length = EEG.pnts * EEG.trials;
else
    data_length = size(EEG.data, 2) * size(EEG.data, 3);
end

if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'ica') && isfield(EEG.ALSUTRECHT.ica, 'icachansind')
    max_avail_chans = length(EEG.ALSUTRECHT.ica.icachansind);
elseif isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    max_avail_chans = length(EEG.icachansind);
else
    max_avail_chans = EEG.nbchan;
end

if isfield(EEG, 'ref') && (strcmpi(EEG.ref, 'averef') || strcmpi(EEG.ref, 'average'))
    max_avail_chans = max_avail_chans - 1;
end

num_ica_max = round(sqrt(data_length / 20));

[~, b] = min(abs(num_ica_candidates - num_ica_max));
selected_n = num_ica_candidates(b);

num_ica = min(selected_n, max_avail_chans);

k_achieved = data_length / (num_ica^2);
dur_min    = (data_length / EEG.srate) / 60;

fprintf('Recording Duration:              %.1f min (%d pts @ %d Hz)\n', dur_min, data_length, EEG.srate);
fprintf('Recommended max{#PCs} (k=20):    %d\n', num_ica_max);
fprintf('Selected #PCs for decomposition: %d (Actual k = %.1f)\n', num_ica, k_achieved);

end

% =========================================================================
% Helper: Variance Tracking
% =========================================================================
function [vaf_compvar, var_global, var_peak_chan, peak_chan_idx] = calc_ic_variance(EEG)

W_inv  = double(EEG.icawinv);
ch_idx = EEG.icachansind;
if isempty(ch_idx); ch_idx = 1:size(W_inv, 1); end

data_2d = double(reshape(EEG.data(ch_idx, :, :), length(ch_idx), []));

if isempty(EEG.icaact)
    icaact = (EEG.icaweights * EEG.icasphere) * data_2d;
else
    icaact = double(reshape(EEG.icaact, size(W_inv, 2), []));
end

num_ics = size(W_inv, 2);
var_act = var(icaact, 0, 2)';

total_raw_var = sum(var(data_2d, 0, 2));
vaf_compvar   = zeros(num_ics, 1);

for i = 1:num_ics
    diff_data = data_2d - (W_inv(:, i) * icaact(i, :));
    vaf_compvar(i) = 100 * (1 - (sum(var(diff_data, 0, 2)) / (total_raw_var + eps)));
end

winv_sq           = W_inv.^2;
proj_var_per_chan = winv_sq .* repmat(var_act, size(W_inv, 1), 1);
total_ic_power    = sum(proj_var_per_chan, 1);
var_global        = 100 * (total_ic_power ./ (sum(total_ic_power) + eps));
var_global        = var_global(:);

var_raw_chans = var(data_2d, 0, 2);
[~, peak_chan_idx] = max(abs(W_inv), [], 1);
var_peak_chan = zeros(num_ics, 1);

for i = 1:num_ics
    c = peak_chan_idx(i);
    var_peak_chan(i) = 100 * (proj_var_per_chan(c, i) / (var_raw_chans(c) + eps));
end

end