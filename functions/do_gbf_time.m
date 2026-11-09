function EEG = do_gbf_time(EEG, cfg)
% DO_GBF_TIME Geometry-Aware Subspace Projection and Time-Aware Burst Suppression
%
% METHOD OVERVIEW:
% This function cleans broadband electromyographic (EMG) contamination
% (both tonic tension and transient clenches/twitches) from continuous 128-channel
% EEG data without distorting low-frequency neural rhythms or creating spectral
% cliffs. It operates via a two-stage hybrid framework:
%
% 1. Geometric Basis Function (GBF) Cortical Subspace Projection:
%    - Extracts smooth spatial eigenmodes of the cortical manifold by solving
%      the generalized eigenvalue problem of the discrete Laplace-Beltrami
%      operator (cotangent stiffness S and lumped mass M) across subject-specific,
%      age-matched cortical meshes (left and right hemispheres treated independently).
%    - Projects oriented, surface-normal dipole leadfields onto these spatial
%      harmonics and computes a Maximum A Posteriori (MAP) inverse operator.
%    - A logarithmic prior penalty, Sigma_inv = diag(-beta / log(lambda)),
%      penalises high-spatial-frequency cortical modes while allowing low,
%      spatially smooth dipolar fields to pass unhindered.
%    - Sensor data above 15 Hz are gradually guided into this cortical subspace
%      P_brain via a cosine transition taper (15 to 35 Hz), eliminating
%      the non-dipolar spatial topographies characteristic of tonic cranial EMG.
%
% 2. Residual-Based Time-Aware Burst Gating:
%    - Detects transient bursts (e.g., masseter/temporalis clenches) by tracking
%      the non-cortical residual R_HF = X_HF - P_brain * X_HF. Because cranial
%      muscles severely violate the leadfield geometry, R_HF spikes drastically
%      during clenches.
%    - Computes a localized, top-K channel metric across time to ensure sensitivity
%      to focal unilateral contractions without global dilution across 128 channels.
%    - Computes a smooth sigmoid attenuation gain calibrated against the quiet-epoch
%      baseline RMS. Crucially, this burst gain is applied exclusively to the
%      high-frequency stream (>20 Hz), leaving slow-wave baselines, blinks,
%      and eye movements untouched to prevent vertical step jumps.
%
% INPUTS:
%   EEG        - EEGLAB structure containing continuous .data, .srate, and subject ID.
%   cfg        - Optional configuration structure:
%                .num_modes    : GBF spatial harmonics per hemisphere (default: 60)
%                .target_ratio : Balance between data term (L^T L) and prior (default: 3.0)
%                .burst_thresh : Detection sensitivity in MAD units for top channels (default: 3.0)
%                .time_smooth  : Temporal integration window for residual envelope (default: 60 ms)
%
% OUTPUTS:
%   EEG        - Cleaned EEGLAB structure preserving 1/f spectral continuity.

%   TARGET_RATIO TUNING GUIDE:
%   Balances the forward model fit against the cortical geometric prior:
%     Tr(L_sq) / Tr(Sigma_inv) = target_ratio
%
%   - Higher target_ratio (e.g. 5 to 10):
%       * Forward data fit dominates over the prior penalty.
%       * Higher spatial-frequency cortical modes are retained with less damping.
%       * Spatial filter P_brain becomes more permissive (effective rank increases;
%         diagonal channel gain rises).
%       * Trade-off: Risk of overfitting to sensor noise or leaking focal
%         peripheral myogenic patterns through complex mode combinations.
%
%   - Lower target_ratio (e.g. 1 to 2):
%       * Geometric penalty dominates over the data fit.
%       * Finer cortical spatial variations are heavily suppressed, preserving
%         only broad, low-frequency Laplace-Beltrami harmonics.
%       * Spatial filter P_brain becomes strictly conservative, constraining
%         sensor data to a low-dimensional cortical subspace.
%       * Trade-off: Risk of over-attenuating genuine focal cortical gamma activity
%         alongside myogenic noise.
%
%   - Default (target_ratio = 3.0):
%       * Maintains forward projection gain dominant by a factor of 3 over the prior.
%       * Typically yields an average diagonal sensor gain of 0.18 to 0.25 on
%         128-channel caps, optimal for cranial EMG suppression.

num_modes     = cfg.gbf.num_modes;
target_ratio  = cfg.gbf.target_ratio;
burst_thresh  = cfg.gbf.burst_thresh;
time_smooth   = cfg.gbf.time_smooth;
f_trans_start = cfg.gbf.freq_trans(1);
f_trans_end   = cfg.gbf.freq_trans(2);

% Freq where EMG contamination starts
f_myo = max(50, f_trans_end);

fprintf('\n================================\n');
fprintf('Time-Aware GBF EMG Subspace Filter\n');
fprintf('================================\n');

% 0. Leadfield Preparation & Subject Matching
% Table1 = readtable('C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-10-13.txt', 'VariableNamingRule', 'preserve');
Table1 = readtable(EEG.ALSUTRECHT.subject.table1, 'VariableNamingRule', 'preserve');
Table1.ID = categorical(Table1.ID);
age = Table1.AGE(Table1.ID == EEG.ALSUTRECHT.subject.id);
age = mean(age);
fprintf('Average age: %.1f\n', age);

age_lf = 20:10:80;
[~, b] = min(abs(age_lf - age));
best_lf = age_lf(b);

file_name = ['leadfield_vrtx_' num2str(best_lf) 'yo.mat'];
file_path = fullfile('C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\template\eeg\BEM\age', file_name);
load(file_path, 'leadfield');

leadfield = orient_leadfield(leadfield);
K_inside = cell2mat(leadfield.leadfield(leadfield.inside).');

indx = find(leadfield.inside, 1, 'first');
num_ori = size(leadfield.leadfield{indx}, 2);
fprintf('Leadfield: %s\n', file_name);
fprintf('Using %d/%d dipoles from the leadfield.\n', size(K_inside, 2), num_ori * sum(leadfield.inside));

% 1. Extract inverse operator using the GBF MAP function
fprintf('1. Calculating GBF spatial filter (%d modes, target_ratio=%g)...\n', num_modes, target_ratio);

% Estimate channel-wise tonic high-frequency noise floor (>50 Hz)
X = double(EEG.data);
X_gamma = highpass(X', f_myo, EEG.srate)';
chan_noise_var = median(X_gamma.^2, 2);
chan_noise_var = chan_noise_var / median(chan_noise_var); % Normalise relative to median scalp channel

[W_inv, W_theta] = do_gbf_map(K_inside, leadfield, num_modes, target_ratio, chan_noise_var);

% 2. Reconstruct the spatial projection matrix (P_brain)
fprintf('2. Constructing sensor-to-sensor projection matrix...\n');
row_norms = vecnorm(K_inside, 2, 2);
K_norm = K_inside ./ row_norms;

% Pure MAP projection (no artificial scaling or clamping)
P_brain = diag(row_norms) * K_norm * W_inv * diag(1 ./ row_norms);

mean_gain = mean(diag(P_brain));
fprintf('Average diagonal channel gain of P_brain: %.2f\n', mean_gain);

% -------------------------------------------------------------------------
% 3. Frequency-Domain Spatial Filtering (Optimized Block-Wise)
% -------------------------------------------------------------------------
fprintf('3. Computing block-wise frequency-dependent subspace projection...\n');

bounds = [];
if isfield(EEG, 'event') && ~isempty(EEG.event)
    for i = 1:length(EEG.event)
        if ischar(EEG.event(i).type) && strcmp(EEG.event(i).type, 'boundary')
            bounds = [bounds, round(EEG.event(i).latency)];
        end
    end
end

seg_starts = [1, bounds];
seg_ends   = [bounds - 1, EEG.pnts];
n_blocks   = length(seg_starts);

X_recon = zeros(size(X), 'single');
last_n_pad_pnts = -1; % Cache tracker for segment length
w_full = [];          % Pre-allocate

for b = 1:n_blocks
    idx_start = seg_starts(b);
    idx_end   = seg_ends(b);
    if idx_end <= idx_start, continue; end

    X_seg = double(X(:, idx_start:idx_end));
    n_seg_pnts = size(X_seg, 2);

    pad_pnts = min(round(2.0 * EEG.srate), n_seg_pnts - 1);
    if pad_pnts < 1
        X_pad = X_seg;
    else
        X_pad = [fliplr(X_seg(:, 1:pad_pnts)), X_seg, fliplr(X_seg(:, end-pad_pnts+1:end))];
    end

    n_pad_pnts = size(X_pad, 2);

    % Only recalculate the spatial taper if the block size changes
    if n_pad_pnts ~= last_n_pad_pnts
        freqs_fft = (0:floor(n_pad_pnts/2)) * (EEG.srate / n_pad_pnts);
        half_len  = numel(freqs_fft);
        w_spatial = zeros(1, half_len);

        in_trans  = (freqs_fft >= f_trans_start) & (freqs_fft <= f_trans_end);
        w_spatial(in_trans) = 0.5 * (1 - cos(pi * (freqs_fft(in_trans) - f_trans_start) / (f_trans_end - f_trans_start)));
        w_spatial(freqs_fft > f_trans_end) = 1.0;

        if mod(n_pad_pnts, 2) == 0
            w_full = [w_spatial, fliplr(w_spatial(2:end-1))];
        else
            w_full = [w_spatial, fliplr(w_spatial(2:end))];
        end
        last_n_pad_pnts = n_pad_pnts; % Update cache
    end

    X_fft       = fft(X_pad, [], 2);
    P_X_fft     = P_brain * X_fft;
    X_clean_fft = (X_fft .* (1 - w_full)) + (P_X_fft .* w_full);

    X_recon_pad = real(ifft(X_clean_fft, [], 2));

    if pad_pnts > 0
        X_recon(:, idx_start:idx_end) = X_recon_pad(:, pad_pnts+1 : end-pad_pnts);
    else
        X_recon(:, idx_start:idx_end) = X_recon_pad;
    end
end

% -------------------------------------------------------------------------
% 4. Global Burst Suppression (Relative Spatial Error Gating)
% -------------------------------------------------------------------------
bounds = [];
if isfield(EEG, 'event') && ~isempty(EEG.event)
    for i = 1:length(EEG.event)
        if ischar(EEG.event(i).type) && strcmp(EEG.event(i).type, 'boundary')
            bounds = [bounds, round(EEG.event(i).latency)];
        end
    end
end

seg_starts = [1, bounds];
seg_ends   = [bounds - 1, EEG.pnts];
n_blocks   = length(seg_starts);

w_samples = max(1, round((time_smooth / 1000) * EEG.srate));
n_chan = size(EEG.data, 1);

% 1. Create a zero-phase Butterworth filter (highly stable, no phase shifts)
[b_myo, a_myo] = butter(4, f_myo / (EEG.srate / 2), 'high');

residual_rms_full = zeros(1, EEG.pnts, 'single');
X_recon_HF_full   = zeros(n_chan, EEG.pnts, 'single');

fprintf('4a. Computing global spatial residual RMS across %d continuous block(s)...\n', n_blocks);

for b = 1:n_blocks
    idx_start = seg_starts(b);
    idx_end   = seg_ends(b);
    if idx_end <= idx_start, continue; end

    X_seg = double(X(:, idx_start:idx_end));
    X_recon_seg = double(X_recon(:, idx_start:idx_end));
    n_seg_pnts = size(X_seg, 2);

    % Pad to prevent edge ringing
    pad_pnts = min(round(1.0 * EEG.srate), n_seg_pnts - 1);
    if pad_pnts > 0
        X_pad = [fliplr(X_seg(:, 1:pad_pnts)), X_seg, fliplr(X_seg(:, end-pad_pnts+1:end))];
        X_recon_pad = [fliplr(X_recon_seg(:, 1:pad_pnts)), X_recon_seg, fliplr(X_recon_seg(:, end-pad_pnts+1:end))];
    else
        X_pad = X_seg;
        X_recon_pad = X_recon_seg;
    end

    % Zero-phase extraction of high frequencies
    X_HF_pad = filtfilt(b_myo, a_myo, X_pad')';
    X_recon_HF_pad = filtfilt(b_myo, a_myo, X_recon_pad')';

    if pad_pnts > 0
        X_HF_seg = X_HF_pad(:, pad_pnts+1 : end-pad_pnts);
        X_recon_HF_seg = X_recon_HF_pad(:, pad_pnts+1 : end-pad_pnts);
    else
        X_HF_seg = X_HF_pad;
        X_recon_HF_seg = X_recon_HF_pad;
    end

    % Residual is strictly the high-frequency energy rejected by the spatial filter
    R_HF_seg = X_HF_seg - X_recon_HF_seg;

    % Global RMS combines all 128 channels into a single powerful 1D metric
    R_rms_seg = sqrt(movmean(sum(R_HF_seg.^2, 1), w_samples));

    residual_rms_full(idx_start:idx_end) = R_rms_seg;
    X_recon_HF_full(:, idx_start:idx_end) = X_recon_HF_seg;
end

fprintf('4b. Estimating robust global baseline statistics...\n');
% Calculate baseline ONLY from the quietest 50% of the data.
% This prevents dense jaw clenches from inflating the baseline and masking themselves.
sorted_rms = sort(double(residual_rms_full));
quiet_idx  = 1 : max(1, round(0.50 * length(sorted_rms)));
res_med = median(sorted_rms(quiet_idx));
res_mad = median(abs(sorted_rms(quiet_idx) - res_med)) * 1.4826;

res_z = (double(residual_rms_full) - res_med) ./ max(res_mad, eps);

% Calculate clean cortical RMS for the attenuation scale factor
clean_mask = res_z < 2.0;
if sum(clean_mask) < 0.05 * numel(clean_mask), clean_mask = true(1, numel(clean_mask)); end
base_rms_global = sqrt(mean(double(X_recon_HF_full(:, clean_mask)).^2, 2));

fprintf('4c. Applying soft-gating suppression strictly within blocks...\n');
data_clean   = zeros(size(X_recon), 'single');
burst_active = false(1, EEG.pnts);
thresh       = burst_thresh;

for b = 1:n_blocks
    idx_start = seg_starts(b);
    idx_end   = seg_ends(b);
    if idx_end <= idx_start, continue; end

    res_z_seg = res_z(idx_start:idx_end);
    HF_seg    = double(X_recon_HF_full(:, idx_start:idx_end));
    X_recon_seg = double(X_recon(:, idx_start:idx_end));

    % Perfect complementary low-frequency split. Guaranteed to cause zero phase distortion.
    LF_seg = X_recon_seg - HF_seg;

    % Continuous sigmoid weighting directly from the Z-score
    steepness = 1.0;
    raw_weight = 1 ./ (1 + exp(-(res_z_seg - thresh) / steepness));

    smooth_pnts = round(0.100 * EEG.srate);
    burst_weight = movmean(raw_weight, smooth_pnts);

    burst_active_seg = burst_weight > 0.5;

    inst_rms_seg = sqrt(movmean(HF_seg.^2, w_samples, 2));
    scale_factor = min(1.0, (base_rms_global * 1.8) ./ max(inst_rms_seg, eps));

    % Apply smooth attenuation only where burst_weight is active
    sample_gain = (1 - burst_weight) + (burst_weight .* scale_factor);

    HF_suppressed_seg = HF_seg .* sample_gain;
    data_clean(:, idx_start:idx_end) = LF_seg + HF_suppressed_seg;
    burst_active(idx_start:idx_end)  = burst_active_seg;
end

% Removed the premature EEG.data = data_clean; here
fprintf('   Flagged massive burst samples: %.2f%% of total recording.\n', mean(burst_active) * 100);

% -------------------------------------------------------------------------
% 5. Check
% -------------------------------------------------------------------------
% Check 1
EEG_NEW = EEG;
EEG_NEW.data = data_clean;

% cut_off = 1;
% EEG_PLOT     = do_filtering_method(EEG, cut_off, cfg);
% EEG_NEW_PLOT = do_filtering_method(EEG_NEW, cut_off, cfg);
% vis_artifacts(EEG_NEW_PLOT, EEG_PLOT);

% Plot 1
check_power_cleaning(EEG, EEG_NEW, [], 'gbf', cfg);

% Plot 2
cfg_diag.freq_range  = [1 100];
cfg_diag.n_fft       = 2;
cfg_diag.plot_diff   = true;
cfg_diag.bands.Theta = [4 8];
cfg_diag.bands.Alpha = [8 13];
cfg_diag.bands.Beta  = [13 30];
cfg_diag.bands.Gamma = [30 100];
plot_gbf_diagnostics(EEG, EEG_NEW, cfg_diag);

% Plot 3
plot_gbf_diagnostics_2(EEG, EEG_NEW, W_theta);

% -------------------------------------------------------------------------
% 6. Store
% -------------------------------------------------------------------------
EEG.data = data_clean;
EEG.ALSUTRECHT.gbf.burst_mask   = burst_active;
EEG.ALSUTRECHT.gbf.burst_weight = burst_weight;

fprintf('GBF cleaning complete.\n');
end

function [W, W_theta, A_inside, lambda] = do_gbf_map(K_inside, leadfield, num_modes, target_ratio, chan_noise_var)
assert(strcmpi(leadfield.unit, 'mm'));
vertices = leadfield.pos;
faces    = leadfield.tri;
inside   = leadfield.inside;

row_norms = vecnorm(K_inside, 2, 2);
K_norm = K_inside ./ row_norms;

n_total = size(vertices, 1);
n_half = n_total / 2;

if mod(n_total, 2) ~= 0
    error('Vertex count is not even; cannot split into equal halves.');
end

idx1 = 1:n_half;
idx2 = (n_half + 1):n_total;

mean_x1 = mean(vertices(idx1, 1));
mean_x2 = mean(vertices(idx2, 1));

if mean_x1 < mean_x2
    fprintf('Detected L/R indices: 1:%d are left, %d:%d are right.\n', n_half, n_half+1, n_total);
else
    fprintf('Detected L/R indices: 1:%d are right, %d:%d are left.\n', n_half, n_half+1, n_total);
end

faces1_mask = all(ismember(faces, idx1), 2);
faces2_mask = all(ismember(faces, idx2), 2);

faces1 = faces(faces1_mask, :);
faces2 = faces(faces2_mask, :) - n_half;

% Min eig value
sigma = -1e-8;

fprintf('Processing left hemisphere...\n');
[S1, M1] = construct_cotangent_laplacian(vertices(idx1, :), faces1);
[A1, D1] = eigs(S1, M1, num_modes, sigma);

fprintf('Processing right hemisphere...\n');
[S2, M2] = construct_cotangent_laplacian(vertices(idx2, :), faces2);
[A2, D2] = eigs(S2, M2, num_modes, sigma);

lambda1 = abs(diag(D1));
lambda2 = abs(diag(D2));

lambda1(1) = lambda1(2);
lambda2(1) = lambda2(2);

[lambda1, sort_idx1] = sort(lambda1, 'ascend');
A1 = A1(:, sort_idx1);

[lambda2, sort_idx2] = sort(lambda2, 'ascend');
A2 = A2(:, sort_idx2);

lambda1 = lambda1 / (max(lambda1) * 1.05);
lambda2 = lambda2 / (max(lambda2) * 1.05);
lambda  = [lambda1; lambda2];

assert(all(lambda < 1 & lambda > 0), 'Eigenvalues must be strictly between 0 and 1');
fprintf('Normalised eigenvalue spectrum: min = %.4e, max = %.4e\n', min(lambda), max(lambda));

A_full = zeros(n_total, num_modes * 2);
A_full(idx1, 1:num_modes) = A1;
A_full(idx2, num_modes + 1:end) = A2;
A_inside = A_full(inside, :);

% % Plot modes
% plot_cortical_mode(vertices, faces, A_full, 30);

% L = K_norm * A_inside;
% L_sq = L' * L;
%
% P_shape = diag(-1 ./ log(lambda));
% beta = (trace(L_sq) / trace(P_shape)) / target_ratio;
% Sigma_inv = diag(-beta ./ log(lambda));
%
% fprintf('Computing the MAP estimator...\n');
% fprintf('Ratio L_sq to Prior: %.2f\n', trace(L_sq) / trace(Sigma_inv));
% W_theta = (L_sq + Sigma_inv) \ L';
%
% W = A_inside * W_theta;

% Weight leadfield by inverse sensor standard deviation
W_e = diag(1 ./ sqrt(chan_noise_var));
L = W_e * K_norm * A_inside;
L_sq = L' * L;

P_shape = diag(-1 ./ log(lambda));
beta = (trace(L_sq) / trace(P_shape)) / target_ratio;
Sigma_inv = diag(-beta ./ log(lambda));

fprintf('Computing the MAP estimator...\n');
fprintf('Ratio L_sq to Prior: %.2f\n', trace(L_sq) / trace(Sigma_inv));

% MAP estimation with sensor weighting
W_theta = (L_sq + Sigma_inv) \ (L' * W_e);
W = A_inside * W_theta;
end

function [S, M] = construct_cotangent_laplacian(vertices, faces)
num_vertices = size(vertices, 1);
num_faces = size(faces, 1);

I = zeros(num_faces * 6, 1);
J = zeros(num_faces * 6, 1);
V = zeros(num_faces * 6, 1);

v1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
v2 = vertices(faces(:, 3), :) - vertices(faces(:, 2), :);
v3 = vertices(faces(:, 1), :) - vertices(faces(:, 3), :);

cross_prod = cross(v1, -v3, 2);
triangle_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));

cot1 = dot(-v3, v1, 2) ./ sqrt(sum(cross(-v3, v1, 2).^2, 2));
cot2 = dot(-v1, v2, 2) ./ sqrt(sum(cross(-v1, v2, 2).^2, 2));
cot3 = dot(-v2, v3, 2) ./ sqrt(sum(cross(-v2, v3, 2).^2, 2));

idx = 1;
for i = 1:3
    i1 = i;
    i2 = mod(i, 3) + 1;
    i3 = mod(i + 1, 3) + 1;

    edges1 = faces(:, i2);
    edges2 = faces(:, i3);

    if i == 1
        cot_weight = cot1;
    elseif i == 2
        cot_weight = cot2;
    else
        cot_weight = cot3;
    end

    I(idx:idx+num_faces-1) = edges1;
    J(idx:idx+num_faces-1) = edges2;
    V(idx:idx+num_faces-1) = -0.5 * cot_weight;
    idx = idx + num_faces;

    I(idx:idx+num_faces-1) = edges2;
    J(idx:idx+num_faces-1) = edges1;
    V(idx:idx+num_faces-1) = -0.5 * cot_weight;
    idx = idx + num_faces;
end

S = sparse(I(1:idx-1), J(1:idx-1), V(1:idx-1), num_vertices, num_vertices);
S = S - spdiags(sum(S, 2), 0, num_vertices, num_vertices);

vertex_areas = accumarray(faces(:), repmat(triangle_areas / 3, 3, 1), [num_vertices, 1]);
M = spdiags(vertex_areas, 0, num_vertices, num_vertices);
end

function lf_out = orient_leadfield(lf_in)
lf_out = lf_in;
inside_idx = find(lf_in.inside);

for i = 1:numel(inside_idx)
    pos_idx = inside_idx(i);
    L = lf_in.leadfield{pos_idx};
    u = lf_in.ori{pos_idx};

    if isempty(L) || isempty(u)
        continue;
    end

    u = u(:);
    norm_u = norm(u);
    if norm_u > 0
        u = u / norm_u;
    end
    lf_out.leadfield{pos_idx} = L * u;
end

if isfield(lf_out, 'leadfielddimord')
    lf_out.leadfielddimord = '{pos}_chan';
end
end

function plot_gbf_diagnostics(EEG_raw, EEG_clean, cfg_plot)
if nargin < 3, cfg_plot = struct(); end
if ~isfield(cfg_plot, 'total_range'), cfg_plot.total_range = [2 80]; end
if ~isfield(cfg_plot, 'freq_range'),  cfg_plot.freq_range  = [1 85]; end
if ~isfield(cfg_plot, 'n_fft'),       cfg_plot.n_fft       = 2; end
if ~isfield(cfg_plot, 'plot_diff'),   cfg_plot.plot_diff   = true; end

if ~isfield(cfg_plot, 'bands')
    cfg_plot.bands.Theta = [4 8];
    cfg_plot.bands.Alpha = [8 13];
    cfg_plot.bands.Beta  = [13 30];
    cfg_plot.bands.Gamma = [30 80];
end

fs = EEG_raw.srate;
nperseg  = round(cfg_plot.n_fft * fs);
noverlap = round(nperseg * 0.5);

X_raw   = double(EEG_raw.data);
X_clean = double(EEG_clean.data);

fprintf('Computing power spectral densities...\n');
[psd_ch_raw, freqs] = pwelch(X_raw', nperseg, noverlap, nperseg, fs);
[psd_ch_clean, ~]   = pwelch(X_clean', nperseg, noverlap, nperseg, fs);

psd_ch_raw   = psd_ch_raw';
psd_ch_clean = psd_ch_clean';

mean_psd_raw   = 10 * log10(mean(psd_ch_raw, 1));
mean_psd_clean = 10 * log10(mean(psd_ch_clean, 1));

f_idx = (freqs >= cfg_plot.freq_range(1)) & (freqs <= cfg_plot.freq_range(2));
freq_sub = freqs(f_idx);

% figure('Name', 'GBF Diagnostics: Power Spectra', 'Color', 'w', 'Position', [100 100 950 420]);
fh = figure('Color', 'w');
freq_sub = log10(freq_sub);

subplot(1, 2, 1);
plot(freq_sub, mean_psd_raw(f_idx), 'Color', [0.85 0.2 0.2], 'LineWidth', 1.8); hold on;
plot(freq_sub, mean_psd_clean(f_idx), 'Color', [0.15 0.45 0.85], 'LineWidth', 1.8);
grid on; box on;
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Power (10*log_{10} \muV^2/Hz)', 'FontSize', 11);
title('Global Scalp Spectrum (Mean All Chans)', 'FontSize', 12);
legend({'Raw', 'Cleaned'}, 'Location', 'northeast');

subplot(1, 2, 2);
psd_diff_db = mean_psd_raw(f_idx) - mean_psd_clean(f_idx);
area(freq_sub, psd_diff_db, 'FaceColor', [0.85 0.3 0.3], 'FaceAlpha', 0.35, 'EdgeColor', [0.7 0.1 0.1]);
grid on; box on;
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Suppression (\Delta dB)', 'FontSize', 11);
title('Suppressed Power (Raw - Cleaned)', 'FontSize', 12);

% Save
save_figure(fh, EEG_raw.ALSUTRECHT.subject.figures, [EEG_raw.ALSUTRECHT.subject.id '_gbf_1'], [25 10]);

total_idx = (freqs >= cfg_plot.total_range(1)) & (freqs <= cfg_plot.total_range(2));
total_pwr_raw   = sum(psd_ch_raw(:, total_idx), 2);
total_pwr_clean = sum(psd_ch_clean(:, total_idx), 2);

band_names = fieldnames(cfg_plot.bands);
n_bands = numel(band_names);
n_rows  = 2 + double(cfg_plot.plot_diff);

cmap_pos  = brewermap([], 'Reds');
cmap_diff = brewermap([], '*RdBu');

% fh = figure('Name', 'GBF Diagnostics: Normalised Band Power Topographies', ...
%     'Color', 'w', 'Position', [120 120 250 * n_bands 230 * n_rows]);
fh = figure('Color', 'w');

for b = 1:n_bands
    b_name = band_names{b};
    b_lims = cfg_plot.bands.(b_name);
    idx_band = (freqs >= b_lims(1)) & (freqs <= b_lims(2));

    % Absolute band power (sum across frequency bins)
    abs_pwr_raw   = sum(psd_ch_raw(:, idx_band), 2);
    abs_pwr_clean = sum(psd_ch_clean(:, idx_band), 2);

    % Relative band power (%)
    rel_pwr_raw   = (abs_pwr_raw ./ total_pwr_raw) * 100;
    rel_pwr_clean = (abs_pwr_clean ./ total_pwr_clean) * 100;

    shared_c_min = 0;
    shared_c_max = max([rel_pwr_raw; rel_pwr_clean]) * 1.05;

    % Row 1: Raw Relative Power
    ax1 = subplot(n_rows, n_bands, b);
    mytopoplot(rel_pwr_raw, [], [], ax1, [shared_c_min shared_c_max]);
    title(sprintf('%s [%d to %d Hz]\nRaw (%%)', b_name, b_lims(1), b_lims(2)), 'FontSize', 11);
    colormap(ax1, cmap_pos);
    c1 = colorbar;
    ylabel(c1, '%', 'FontSize', 9);

    % Row 2: Cleaned Relative Power
    ax2 = subplot(n_rows, n_bands, n_bands + b);
    mytopoplot(rel_pwr_clean, [], [], ax2, [shared_c_min shared_c_max]);
    title(sprintf('%s\nCleaned (%%)', b_name), 'FontSize', 11);
    colormap(ax2, cmap_pos);
    c2 = colorbar;
    ylabel(c2, '%', 'FontSize', 9);

    % Row 3: Absolute Power Difference (Cleaned - Raw in \muV^2)
    if cfg_plot.plot_diff
        ax3 = subplot(n_rows, n_bands, 2 * n_bands + b);
        diff_abs = abs_pwr_clean - abs_pwr_raw;

        max_abs = max(abs(diff_abs));
        if max_abs == 0, max_abs = eps; end

        mytopoplot(diff_abs, [], [], ax3, [-max_abs max_abs]);
        title(sprintf('%s\n\\Delta Absolute Power', b_name), 'FontSize', 11);
        colormap(ax3, cmap_diff);
        c3 = colorbar;
        ylabel(c3, '\muV^2', 'FontSize', 9);
    end
end

% Save
save_figure(fh, EEG_raw.ALSUTRECHT.subject.figures, [EEG_raw.ALSUTRECHT.subject.id '_gbf_2'], [25 15]);

fprintf('Diagnostic figures successfully generated.\n');
end

function plot_cortical_mode(vertices, faces, mode_data, mode_number)
% PLOT_CORTICAL_MODE Visualises a specific geometric basis function
% on the 3D cortical surface mesh for debugging.
%
% Inputs:
%   vertices    : Spatial coordinates of the cortical surface (Sources by 3)
%   faces       : Triangulation connectivity of the surface (Triangles by 3)
%   mode_data   : The spatial pattern of the mode to plot (Sources by modes)
%   mode_number : Integer identifying the mode number for the title

mode_data = mode_data(:, mode_number);

% Create a new figure with a white background
figure('Name', sprintf('Cortical Mode %d', mode_number), 'Color', 'w');

% Plot the 3D mesh using trisurf
% 'FaceColor', 'interp' smoothly interpolates the colours across the triangles
trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
    mode_data, 'EdgeColor', 'none', 'FaceColor', 'interp');

% Optimise the camera view and axis proportions
axis equal;
axis off;
view([-90, 0]); % Default to a lateral view; adjust angles as needed

% Add lighting and shading to enhance the 3D curvature of the brain folds
camlight('headlight');
lighting gouraud;
material dull;

% Add a colourbar and set the colour map
% 'jet' or standard diverging colour maps work well for oscillating modes
colorbar;
colormap(gca, 'jet');

% Centre the colour limits around zero
% This ensures that negative values and positive values are symmetrical
max_val = max(abs(mode_data));
if max_val > 0
    clim([-max_val, max_val]); % Use caxis([-max_val, max_val]) for older MATLAB versions
end

% Add an informative title
title(sprintf('Geometric Basis Function: Mode %d', mode_number), 'FontSize', 14);
end


function plot_mask(EEG, burst_mask)
% Convert logical mask vector into [start_sample, end_sample] windows
diff_mask = diff([0, burst_mask(:)', 0]);
onsets  = find(diff_mask == 1);
offsets = find(diff_mask == -1) - 1;

n_chans = size(EEG.data, 1);
n_bursts = numel(onsets);

% Construct winrej matrix: [start end R G B ch1 ch2 ... chN]
winrej = zeros(n_bursts, 5 + n_chans);
winrej(:, 1) = onsets(:);
winrej(:, 2) = offsets(:);
winrej(:, 3:5) = repmat([0.9 0.4 0.4], n_bursts, 1); % Soft red colour
winrej(:, 6:end) = 1; % Apply the highlight across all channels

% Launch interactive scrolling viewer
eegplot(EEG.data, 'srate', EEG.srate, ...
    'title', sprintf('Flagged EMG Bursts (Total: %d)', n_bursts), ...
    'eloc_file', EEG.chanlocs, ...
    'winrej', winrej, ...
    'winlength', 10, ...
    'color', 'off');
end

function plot_gbf_diagnostics_2(EEG_raw, EEG_clean, W_theta)
% PLOT_GBF_REMOVED_DIAGNOSTICS Diagnostic inspection of rejected GBF variance.
%
% Inputs:
%   EEG_raw   : Raw EEGLAB struct or 2D matrix [128 x n_samples]
%   EEG_clean : Cleaned EEGLAB struct or 2D matrix [128 x n_samples]
%   W_theta   : GBF MAP inverse projection matrix [n_modes x 128]

% 1. Extract 2D sensor matrices
if isstruct(EEG_raw)
    X_raw = double(EEG_raw.data(1:128, :));
else
    X_raw = double(EEG_raw(1:128, :));
end

if isstruct(EEG_clean)
    X_clean = double(EEG_clean.data(1:128, :));
else
    X_clean = double(EEG_clean(1:128, :));
end

% Ensure 2D continuous array if data was epoched [chan x pnts x trials]
if ndims(X_raw) == 3
    X_raw   = reshape(X_raw, 128, []);
    X_clean = reshape(X_clean, 128, []);
end

% Isolate rejected signal
X_removed = X_raw - X_clean;
n_pnts    = size(X_removed, 2);

% =========================================================================
% Part A: Mode Power Projection (Laplace-Beltrami Harmonics)
% =========================================================================
% Theta_removed: [n_modes x n_pnts]
Theta_removed = double(W_theta) * X_removed;

% Mean power (variance) per spatial mode
mode_pwr = mean(Theta_removed.^2, 2);
n_modes  = length(mode_pwr);

% =========================================================================
% Part B: Spatial PCA on the Sensor-Space Residual
% =========================================================================
% Spatial sample covariance
C_rem = (X_removed * X_removed') / (n_pnts - 1);

[V, D] = eig(C_rem);
[eig_vals, sort_idx] = sort(diag(D), 'descend');
V_sorted = V(:, sort_idx);

% Matrix rank and variance thresholds
tol          = length(eig_vals) * eps(eig_vals(1));
mat_rank     = sum(eig_vals > tol);
var_exp      = (eig_vals / sum(eig_vals)) * 100;
cum_var      = cumsum(var_exp);
opt_cutoff   = find(cum_var >= 90, 1, 'first');
if isempty(opt_cutoff), opt_cutoff = mat_rank; end

fprintf('\n--- GBF Rejection Diagnostics ---\n');
fprintf('Sensor-space rank of removed signal: %d / 128\n', mat_rank);
fprintf('Recommended PCA component cutoff (90%% variance): %d\n', opt_cutoff);
fprintf('First 10 PCs explain: %.1f%% of total removed variance\n\n', sum(var_exp(1:10)));

% =========================================================================
% Part C: Visualisation (3-Row Layout)
% =========================================================================
figure('Name', 'GBF Rejected Subspace: Modes & PCA Topographies', ...
    'Color', 'w', 'Position', [100, 80, 1200, 780]);

% --- Row 1: Mode Power Distribution Bar Plot ---
subplot(3, 5, 1:5);
b_plot = bar(1:n_modes, mode_pwr, 'FaceColor', [0.25 0.45 0.75], 'EdgeColor', 'none');
hold on; grid on; box on;

% Highlight low-order lobar modes (1 to 20) in a distinct colour
if n_modes >= 20
    b_low = bar(1:20, mode_pwr(1:20), 'FaceColor', [0.85 0.35 0.25], 'EdgeColor', 'none');
    legend([b_low, b_plot], {'Low Modes (1-20, Lobar Scale)', 'High Modes (>20, Gyral Scale)'}, ...
        'Location', 'northeast', 'FontSize', 9);
end

xlabel('Laplace-Beltrami Mode Index (k)', 'FontSize', 11);
ylabel('Mean Power (\theta_k^2)', 'FontSize', 11);
title(sprintf('Removed Power Projected into Mode Space (K = %d Modes)', n_modes), 'FontSize', 12);
xlim([0.5, n_modes + 0.5]);

% --- Rows 2 & 3: Removed PCA Topographies (1 to 10) ---
cmap_div = brewermap([], '*RdBu');

for i = 1:10
    % Subplot index: Row 2 uses slots 6 to 10; Row 3 uses slots 11 to 15
    ax = subplot(3, 5, 5 + i);

    topo = V_sorted(:, i);

    % Enforce positive polarity on the dominant peak for consistency
    [~, max_idx] = max(abs(topo));
    if topo(max_idx) < 0
        topo = -topo;
    end

    c_lim = max(abs(topo));
    if c_lim == 0, c_lim = 1; end

    % topoplot(topo, chanlocs, 'maplimits', [-c_lim, c_lim], 'electrodes', 'off');
    mytopoplot(topo, [], [], ax, [-c_lim, c_lim]);
    colormap(ax, cmap_div);

    % Flag whether this component falls within the 90% variance threshold
    if i <= opt_cutoff
        flag_str = '[<=90%]';
    else
        flag_str = '';
    end

    title(sprintf('Removed PC %d %s\n%.1f%% (Cum: %.1f%%)', ...
        i, flag_str, var_exp(i), cum_var(i)), ...
        'FontSize', 9.5, 'FontWeight', 'bold', 'Interpreter', 'none');
end
end