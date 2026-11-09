function sim_results = simulate_gbf_validation(EEG, cfg)
% SIMULATE_GBF_VALIDATION Ground-truth synthetic cortical dipole injection
% to benchmark GBF subspace filtering against standard temporal filtering.

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'target_freq'),  cfg.target_freq   = 45;   end % 40 Hz gamma
if ~isfield(cfg, 'target_snr_db'),cfg.target_snr_db = -6;  end % -6 dB (submerged in noise)

fprintf('\n===============================================\n');
fprintf('Running Ground-Truth Cortical Signal Simulation\n');
fprintf('===============================================\n');

% -------------------------------------------------------------------------
% 1. Load Subject Leadfield & Geodesic Gaussian Patch Construction
% -------------------------------------------------------------------------
file_name = 'ALS31269_3_leadfield_vrtx.mat';
file_path = fullfile('C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\individual\backup\leadfields', file_name);
load(file_path, 'leadfield');

leadfield = orient_leadfield(leadfield);
K_inside = cell2mat(leadfield.leadfield(leadfield.inside).');

vertices = leadfield.pos;
faces    = leadfield.tri;
inside   = leadfield.inside;

% Target seed on left sensorimotor cortex (approx: X=-35, Y=-15, Z=55)
target_coord = [-35, -15, 55];
euc_dists = vecnorm(vertices - target_coord, 2, 2);

inside_indices = find(inside);
[~, nearest_idx] = min(euc_dists(inside_indices));
seed_vertex = inside_indices(nearest_idx); % Index into the full vertex set

% Build symmetric mesh graph with Euclidean edge lengths
edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
edge_lens = vecnorm(vertices(edges(:, 1), :) - vertices(edges(:, 2), :), 2, 2);

% Construct sparse graph using undirected connectivity
G = graph(edges(:, 1), edges(:, 2), edge_lens, size(vertices, 1));

% Compute exact geodesic distances across the cortical manifold from the seed
geo_dists_full = distances(G, seed_vertex); % Column vector (all vertices)
geo_dists = geo_dists_full(inside);         % Subselect only valid dipole positions

% Define geodesic patch width (12 to 15 mm geodesic spread)
patch_sigma = 15.0;

% Evaluate geodesic Gaussian weights as a column vector
patch_weights = exp(-(geo_dists(:).^2) / (2 * patch_sigma^2));
patch_weights(geo_dists(:) > (3 * patch_sigma)) = 0;

K_patch  = double(K_inside) * patch_weights;
k_dipole = K_patch / norm(K_patch);

% Visualise the exact geodesic patch used in the simulation
plot_cortical_patch(leadfield, seed_vertex, geo_dists_full, patch_sigma);

% figure; hold on;
% indx = find(leadfield.inside);
% ft_plot_mesh(leadfield.pos(indx, :), 'vertexcolor', 'k');
% ft_plot_mesh(leadfield.pos(indx(target_vertex_idx), :), 'vertexcolor', 'r');

% -------------------------------------------------------------------------
% 2. Synthesize Ground-Truth Neural Waveform (40 Hz Amplitude Modulated)
% -------------------------------------------------------------------------
fs = EEG.srate;
n_pnts = EEG.pnts;
t = (0:(n_pnts - 1)) / fs;

% 40 Hz burst with 2 Hz physiological envelope
f_sig = cfg.target_freq;
s_source = sin(2 * pi * f_sig * t) .* (0.5 * (1 + sin(2 * pi * 2 * t)));

% Project dipole to scalp sensors
X_signal_clean = k_dipole * s_source; % 128 x n_pnts

% -------------------------------------------------------------------------
% 3. Mix with Genuine Scalp EMG at Calibrated SNR
% -------------------------------------------------------------------------
EEG = do_filtering_method(EEG, 1, cfg);
X_noise = double(EEG.data);

pwr_noise  = mean(var(X_noise, 0, 2));
pwr_signal = mean(var(X_signal_clean, 0, 2));

% Scale signal to match desired SNR (in dB)
target_pwr = pwr_noise * (10^(cfg.target_snr_db / 10));
scale_factor = sqrt(target_pwr / pwr_signal);
X_signal_scaled = X_signal_clean * scale_factor;

% Create synthetic contaminated dataset
EEG_synth = EEG;
EEG_synth.data = single(X_noise + X_signal_scaled);

fprintf('Injected %d Hz dipole at vertex %d (Left Sensorimotor, SNR = %.1f dB)\n', ...
    f_sig, nearest_idx, cfg.target_snr_db);

% -------------------------------------------------------------------------
% 4. Run Filter Competitors
% -------------------------------------------------------------------------
% A. Low-pass filter (30 Hz cutoff, 4th order zero-phase Butterworth)
fprintf('Running standard temporal low-pass filter (30 Hz cutoff)...\n');
[b_lp, a_lp] = butter(4, 30 / (fs / 2), 'low');
X_lp = filtfilt(b_lp, a_lp, double(EEG_synth.data)')';

% B. Time-Aware GBF Subspace Filter
fprintf('Running Time-Aware GBF Filter...\n');
cfg_gbf = cfg;
EEG_gbf = do_gbf_time(EEG_synth, cfg_gbf);
X_gbf   = double(EEG_gbf.data);

% -------------------------------------------------------------------------
% 5. Quantitative Benchmarking & SNR Metrics
% -------------------------------------------------------------------------
% Choose central electrode closest to the dipole projection peak
[~, best_chan] = max(abs(k_dipole));
chan_label = EEG.chanlocs(best_chan).labels;

% Target waveform at this sensor
true_sig = X_signal_scaled(best_chan, :);

% Compute Waveform Correlations (Phase and Envelope preservation)
corr_mix = corr(true_sig', EEG_synth.data(best_chan, :)');
corr_lp  = corr(true_sig', X_lp(best_chan, :)');
corr_gbf = corr(true_sig', X_gbf(best_chan, :)');

% Compute Relative Root Mean Square Error (RRMSE)
rrmse = @(est, ref) norm(est - ref) / norm(ref);
rrmse_mix = rrmse(double(EEG_synth.data(best_chan, :)), true_sig);
rrmse_lp  = rrmse(X_lp(best_chan, :), true_sig);
rrmse_gbf = rrmse(X_gbf(best_chan, :), true_sig);

% Compute Output SNR in the frequency bin of interest
nfft = 2 * fs;
[pxx_true, f] = pwelch(true_sig, nfft, fs, nfft, fs);
[pxx_mix, ~]  = pwelch(double(EEG_synth.data(best_chan, :)), nfft, fs, nfft, fs);
[pxx_lp, ~]   = pwelch(X_lp(best_chan, :), nfft, fs, nfft, fs);
[pxx_gbf, ~]  = pwelch(X_gbf(best_chan, :), nfft, fs, nfft, fs);

[~, f_idx] = min(abs(f - f_sig));
snr_out = @(pxx) 10 * log10(pxx(f_idx) / mean(pxx([f_idx-4:f_idx-2, f_idx+2:f_idx+4])));

fprintf('\n--- BENCHMARK RESULTS (Sensor: %s) ---\n', chan_label);
fprintf('Method             | Correlation | RRMSE    | %d Hz Peak SNR (dB)\n', f_sig);
fprintf('----------------------------------------------------------\n');
fprintf('Raw Contaminated   |   %6.3f    | %7.2f  |     %6.2f dB\n', corr_mix, rrmse_mix, snr_out(pxx_mix));
fprintf('Temporal Low-Pass  |   %6.3f    | %7.2f  |     %6.2f dB\n', corr_lp,  rrmse_lp,  snr_out(pxx_lp));
fprintf('Time-Aware GBF     |   %6.3f    | %7.2f  |     %6.2f dB\n', corr_gbf, rrmse_gbf, snr_out(pxx_gbf));
fprintf('----------------------------------------------------------\n');

% -------------------------------------------------------------------------
% 6. Plotting Ground-Truth Recovery
% -------------------------------------------------------------------------
figure('Name', 'GBF vs Low-Pass Ground-Truth Simulation', 'Color', 'w', 'Position', [150 150 1100 650]);

t_crop = (t >= 5) & (t <= 7); % Zoom into a 2-second burst segment

subplot(3, 1, 1);
plot(t(t_crop), true_sig(t_crop), 'k', 'LineWidth', 1.8); hold on;
plot(t(t_crop), EEG_synth.data(best_chan, t_crop), 'Color', [0.8 0.2 0.2 0.6], 'LineWidth', 1.0);
legend({'Ground-Truth Cortical Signal', 'Scalp Contaminated Mixture'}, 'Location', 'northeast');
title(sprintf('Synthetic Mixture at Sensor %s (Input SNR: %.1f dB)', chan_label, cfg.target_snr_db), 'FontSize', 11);
xlabel('Time (s)'); ylabel('\muV'); grid on; box on;

subplot(3, 1, 2);
plot(t(t_crop), true_sig(t_crop), 'k', 'LineWidth', 1.8); hold on;
plot(t(t_crop), X_lp(best_chan, t_crop), 'Color', [0.85 0.55 0.1], 'LineWidth', 1.2);
legend({'Ground-Truth Cortical Signal', 'Standard Low-Pass (30 Hz)'}, 'Location', 'northeast');
title(sprintf('Standard Temporal Low-Pass Filter (Corr: %.3f | Signal Erased)', corr_lp), 'FontSize', 11);
xlabel('Time (s)'); ylabel('\muV'); grid on; box on;

subplot(3, 1, 3);
plot(t(t_crop), true_sig(t_crop), 'k', 'LineWidth', 1.8); hold on;
plot(t(t_crop), X_gbf(best_chan, t_crop), 'Color', [0.15 0.5 0.85], 'LineWidth', 1.2);
legend({'Ground-Truth Cortical Signal', 'Time-Aware GBF Subspace'}, 'Location', 'northeast');
title(sprintf('Time-Aware GBF Reconstruction (Corr: %.3f | Phase & Envelope Preserved)', corr_gbf), 'FontSize', 11);
xlabel('Time (s)'); ylabel('\muV'); grid on; box on;

% Store return struct
sim_results = struct();
sim_results.best_chan   = chan_label;
sim_results.corr_lp     = corr_lp;
sim_results.corr_gbf    = corr_gbf;
sim_results.rrmse_lp    = rrmse_lp;
sim_results.rrmse_gbf   = rrmse_gbf;
sim_results.snr_raw     = snr_out(pxx_mix);
sim_results.snr_lp      = snr_out(pxx_lp);
sim_results.snr_gbf     = snr_out(pxx_gbf);
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

function plot_cortical_patch(leadfield, seed_vertex, geo_dists_full, patch_sigma)
% PLOT_CORTICAL_PATCH Visualises the true geodesic cortical patch on the mesh.
%
% Inputs:
%   leadfield      : FieldTrip-style struct containing .pos, .tri
%   seed_vertex    : Vertex index of the patch seed/centre
%   geo_dists_full : Exact geodesic distance vector from graph shortest-path (all vertices)
%   patch_sigma    : Standard deviation of Gaussian drop-off in mm (e.g. 15.0)

if nargin < 4 || isempty(patch_sigma), patch_sigma = 15.0; end

vertices = leadfield.pos;
faces    = leadfield.tri;

% Evaluate geodesic Gaussian weights across all mesh vertices
patch_data = exp(-(geo_dists_full(:).^2) / (2 * patch_sigma^2));
patch_data(geo_dists_full(:) > (3 * patch_sigma)) = 0; % Cut off manifold tails

seed_coord = vertices(seed_vertex, :);

% Render 3D cortical surface
figure('Name', 'Simulated Cortical Patch (Geodesic)', 'Color', 'w', 'Position', [150 150 950 700]);

% Render mesh with interpolated vertex weights
trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
    patch_data, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;

% Mark the exact seed vertex
plot3(seed_coord(1), seed_coord(2), seed_coord(3), 'k.', 'MarkerSize', 24);

% Proportions, view, and lighting
axis equal;
axis off;

% Choose lateral perspective based on seed hemisphere (X coordinate)
if seed_coord(1) < 0
    view([-85, 20]);  % Left lateral view
else
    view([85, 20]);   % Right lateral view
end

camlight('headlight');
lighting gouraud;
material([0.3 0.8 0.1 10 1.0]); % Clean diffuse surface response

% Cortical colormap: neutral anatomical grey to warm activation colours
n_bins = 64;
grey_base = repmat([0.82 0.82 0.84], round(n_bins * 0.25), 1);
hot_ramp  = [linspace(0.82, 0.95, n_bins)', linspace(0.82, 0.20, n_bins)', linspace(0.84, 0.10, n_bins)'];
colormap(gca, [grey_base; hot_ramp]);
clim([0, 1]);

cb = colorbar;
ylabel(cb, 'Normalised Geodesic Weight', 'FontSize', 10);

title(sprintf('Geodesic Cortical Patch (\\sigma = %.1f mm)\nSeed Location: [%.1f, %.1f, %.1f] mm', ...
    patch_sigma, seed_coord(1), seed_coord(2), seed_coord(3)), 'FontSize', 12);
end