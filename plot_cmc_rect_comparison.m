function [cmc_results, fh] = plot_cmc_rect_comparison(DATA, cfg)
% PLOT_CMC_RECT_COMPARISON
% Re-references scalp EEG to Common Average Reference (CAR), extracts
% full task epochs (-2s pre-stim '21' to +3s post-release '22', total [-2.0, 8.0]s),
% and computes 1D stationary CMC (1-5s hold) alongside smoothed Time-Frequency CMC.
%
% Inputs:
%   DATA : EEGLAB dataset structure (sampled at 2048 Hz, continuous or epoched)
%   cfg  : (Optional) Configuration structure:
%          cfg.epoch_win_sec : Full extracted epoch in seconds (default: [-2.0, 8.0])
%          cfg.task_win_sec  : Steady-state hold window for 1D CMC (default: [1.0, 5.0])
%          cfg.event_marker  : Trigger marker for trial onset (default: '21')
%          cfg.method        : 1D spectral method ('multitaper' [default] or 'welch')
%          cfg.nw            : Time-halfbandwidth product for DPSS (default: 3.0)
%          cfg.num_tapers    : Number of Slepian tapers (default: 2*nw - 1 -> 5)
%          cfg.do_tf         : Compute and plot Time-Frequency CMC (default: true)
%          cfg.tf_n_cycles   : Number of Morlet wavelet cycles (default: 5)
%          cfg.apply_car     : Apply Common Average Reference (default: true)
%          cfg.eeg_chan      : Center EEG channel label (default: 'D19')
%          cfg.emg_chan      : EMG muscle name 'APB' or 'FDI' (default: 'APB')
%          cfg.emg_hp        : EMG high-pass cutoff in Hz (default: 5 Hz)
%          cfg.freq_range    : Frequency display range (default: [1, 48])
%          cfg.win_sec       : Welch sub-window duration in seconds (default: 1.0)
%          cfg.overlap_pct   : Welch sub-window overlap percentage (default: 50)
%          cfg.do_plot       : Render summary plot (default: true)
%
% Outputs:
%   cmc_results : Struct containing frequencies, time vectors, 1D and TF coherence arrays
%   fh          : Figure handle

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'epoch_win_sec'), cfg.epoch_win_sec = [-3.0, 8.0];  end % -2s pre-21 to +3s post-22
if ~isfield(cfg, 'task_win_sec'),  cfg.task_win_sec  = [2.0, 5.0];   end % Steady-state hold
if ~isfield(cfg, 'event_marker'),  cfg.event_marker  = '21';         end
if ~isfield(cfg, 'method'),        cfg.method        = 'multitaper'; end
if ~isfield(cfg, 'nw'),            cfg.nw            = 5.0;          end
if ~isfield(cfg, 'num_tapers'),    cfg.num_tapers    = floor(2 * cfg.nw - 1); end
if ~isfield(cfg, 'do_tf'),         cfg.do_tf         = true;         end
if ~isfield(cfg, 'tf_n_cycles'),   cfg.tf_n_cycles   = 5;            end
if ~isfield(cfg, 'apply_car'),     cfg.apply_car     = true;         end
if ~isfield(cfg, 'eeg_chan'),      cfg.eeg_chan      = 'D19';        end
if ~isfield(cfg, 'emg_chan'),      cfg.emg_chan      = 'FDI';        end % APB, FDI, FPB, EPB
if ~isfield(cfg, 'emg_hp'),        cfg.emg_hp        = 10;           end % Hz
if ~isfield(cfg, 'freq_range'),    cfg.freq_range    = [5, 48];      end % Hz
if ~isfield(cfg, 'win_sec'),       cfg.win_sec       = 1.0;          end
if ~isfield(cfg, 'overlap_pct'),   cfg.overlap_pct   = 50;           end
if ~isfield(cfg, 'do_plot'),       cfg.do_plot       = true;         end

fs = DATA.srate; % 2048 Hz
all_labels = {DATA.chanlocs.labels};

% -------------------------------------------------------------------------
% 2. Assign Electrode Coordinates from Template and Label Types
% -------------------------------------------------------------------------
mat_data = load('biosemi128_eeglab.mat', 'chanlocs');
template_locs = mat_data.chanlocs;

data_labels     = {DATA.chanlocs.labels};
template_labels = {template_locs.labels};

coord_fields = intersect(fieldnames(template_locs), ...
    {'X', 'Y', 'Z', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', ...
    'sph_theta_besa', 'sph_phi_besa', 'type', 'ref'});

matched_eeg = false(1, length(DATA.chanlocs));

for i = 1:length(template_locs)
    idx_match = find(strcmpi(data_labels, template_labels{i}), 1);
    if ~isempty(idx_match)
        matched_eeg(idx_match) = true;
        for f = 1:length(coord_fields)
            fld = coord_fields{f};
            DATA.chanlocs(idx_match).(fld) = template_locs(i).(fld);
        end
        if isempty(DATA.chanlocs(idx_match).type)
            DATA.chanlocs(idx_match).type = 'EEG';
        end
    end
end

for k = find(~matched_eeg)
    DATA.chanlocs(k).type = 'EXT';
end

% -------------------------------------------------------------------------
% 3. Common Average Referencing (CAR) across Scalp Channels
% -------------------------------------------------------------------------
scalp_idx = find(strcmpi({DATA.chanlocs.type}, 'EEG'));
if isempty(scalp_idx), scalp_idx = 1:min(128, length(DATA.chanlocs)); end

if cfg.apply_car
    fprintf('Applying Common Average Reference (CAR) across %d scalp channels...\n', length(scalp_idx));
    car_mean = mean(DATA.data(scalp_idx, :), 1);
    DATA.data(scalp_idx, :) = DATA.data(scalp_idx, :) - car_mean;
end

% -------------------------------------------------------------------------
% Fitler
% -------------------------------------------------------------------------
% DATA = do_filtering_fir(DATA);

% -------------------------------------------------------------------------
% 4. Resolve Center EEG Channel and Local Surface Laplacian
% -------------------------------------------------------------------------
idx_center = find(strcmpi(all_labels, cfg.eeg_chan), 1);
assert(~isempty(idx_center), 'Center EEG channel [%s] not found in DATA.', cfg.eeg_chan);
eeg_desc = sprintf('%s (CAR Monopolar)', cfg.eeg_chan);

% -------------------------------------------------------------------------
% 5. Resolve EMG Channel (Bipolar derivation from monopolar EXG leads)
% -------------------------------------------------------------------------
req_chan = cfg.emg_chan;

switch upper(req_chan)
    case 'APB'
        pairs_to_try = {'E1', 'E2'};
        muscle_tag   = 'APB';
    case {'FDI'}
        pairs_to_try = {'E5', 'E6'};
        muscle_tag   = 'FDI';
    case 'FPB'
        pairs_to_try = {'E7', 'E8'};
        muscle_tag   = 'FPB';
    case 'EPB'
        pairs_to_try = {'E9', 'E10'};
        muscle_tag   = 'EPB';
    otherwise
        error('check_emg:UnknownMuscle', 'Unsupported EMG target [%s].', req_chan);
end


idx_pos = find(strcmpi(all_labels, pairs_to_try{1}) | strcmpi(all_labels, ['EXG' pairs_to_try{1}(2:end)]), 1);
idx_neg = find(strcmpi(all_labels, pairs_to_try{2}) | strcmpi(all_labels, ['EXG' pairs_to_try{2}(2:end)]), 1);

assert(~isempty(idx_pos) && ~isempty(idx_neg), ...
    'Could not resolve monopolar pair [%s, %s] for %s.', pairs_to_try{1}, pairs_to_try{2}, muscle_tag);

idx_emg   = idx_pos;
emg_label = sprintf('%s (%s-%s)', muscle_tag, all_labels{idx_pos}, all_labels{idx_neg});

% Bipolar subtraction (Pos - Neg)
DATA.data(idx_emg, :, :) = DATA.data(idx_pos, :, :) - DATA.data(idx_neg, :, :);

if isfield(DATA, 'chanlocs') && numel(DATA.chanlocs) >= idx_emg
    DATA.chanlocs(idx_emg).labels = emg_label;
end

% -------------------------------------------------------------------------
% 6. Full Epoch Extraction ([-2.0s, +8.0s] relative to '21')
% -------------------------------------------------------------------------
win_start_samp = round(cfg.epoch_win_sec(1) * fs); % -4096 pts @ 2048 Hz
win_end_samp   = round(cfg.epoch_win_sec(2) * fs); % +16384 pts @ 2048 Hz
epoch_pnts     = win_end_samp - win_start_samp + 1;
time_vec       = linspace(cfg.epoch_win_sec(1), cfg.epoch_win_sec(2), epoch_pnts);

eeg_epochs = [];
emg_epochs = [];

if ndims(DATA.data) == 3 && isfield(DATA, 'times')
    time_mask = (DATA.times >= cfg.epoch_win_sec(1) * 1000) & (DATA.times <= cfg.epoch_win_sec(2) * 1000);

    eeg_raw = double(DATA.data(idx_center, time_mask, :));
    emg_raw    = double(DATA.data(idx_emg, time_mask, :));
    eeg_epochs = squeeze(eeg_raw);
    emg_epochs = squeeze(emg_raw);
else
    assert(isfield(DATA, 'event') && ~isempty(DATA.event), 'Event structure required for epoching.');

    ev_types   = {DATA.event.type};
    target_str = string(cfg.event_marker);

    ev_matches = false(1, length(ev_types));
    for ev_i = 1:length(ev_types)
        curr_ev = string(ev_types{ev_i});
        if curr_ev == target_str || strcmpi(curr_ev, "S " + target_str) || strcmpi(curr_ev, "S" + target_str)
            ev_matches(ev_i) = true;
        end
    end

    event_latencies = round([DATA.event(ev_matches).latency]);
    assert(~isempty(event_latencies), 'No events matching marker [%s] found.', target_str);

    n_total_pnts = size(DATA.data, 2);

    for k = 1:length(event_latencies)
        t0 = event_latencies(k);
        idx_epoch = (t0 + win_start_samp) : (t0 + win_end_samp);

        if idx_epoch(1) >= 1 && idx_epoch(end) <= n_total_pnts
            eeg_snip = double(DATA.data(idx_center, idx_epoch));
            emg_snip = double(DATA.data(idx_emg, idx_epoch));

            eeg_epochs(:, end+1) = eeg_snip(:); %#ok<AGROW>
            emg_epochs(:, end+1) = emg_snip(:); %#ok<AGROW>
        end
    end
end

% eeg_epochs = eeg_epochs(:, 3:end);
% emg_epochs = emg_epochs(:, 3:end);

n_valid_trials = size(eeg_epochs, 2);
assert(n_valid_trials >= 3, 'Insufficient valid trials (N=%d) extracted for CMC.', n_valid_trials);

% -------------------------------------------------------------------------
% 7. EMG Conditioning (Unrectified vs. Rectified)
% -------------------------------------------------------------------------
[b_hp, a_hp] = butter(4, cfg.emg_hp / (fs / 2), 'high');
emg_filt = filtfilt(b_hp, a_hp, emg_epochs);

emg_unrect = emg_filt - mean(emg_filt, 1);
emg_rect   = abs(emg_filt);
emg_rect   = emg_rect - mean(emg_rect, 1);

eeg_sig = eeg_epochs - mean(eeg_epochs, 1);

% -------------------------------------------------------------------------
% 8. 1D Spectral Estimation on Steady-State Hold ([+1.0s, +5.0s])
% -------------------------------------------------------------------------
hold_mask   = (time_vec >= cfg.task_win_sec(1)) & (time_vec <= cfg.task_win_sec(2));
eeg_hold    = eeg_sig(hold_mask, :);
emg_u_hold  = emg_unrect(hold_mask, :);
emg_r_hold  = emg_rect(hold_mask, :);
n_hold_pnts = size(eeg_hold, 1);

switch lower(cfg.method)
    case 'multitaper'
        [tapers, ~] = dpss(n_hold_pnts, cfg.nw, cfg.num_tapers);

        n_fft = max(2048, 2^nextpow2(n_hold_pnts));
        n_freq_bins = (n_fft / 2) + 1;
        f_axis = (0:(n_fft / 2))' * (fs / n_fft);

        K = cfg.num_tapers;
        total_L = n_valid_trials * K;

        S_ee        = zeros(n_freq_bins, 1);
        S_mm_unrect = zeros(n_freq_bins, 1);
        S_em_unrect = zeros(n_freq_bins, 1);
        S_mm_rect   = zeros(n_freq_bins, 1);
        S_em_rect   = zeros(n_freq_bins, 1);

        for tr = 1:n_valid_trials
            eeg_tr   = eeg_hold(:, tr);
            emg_u_tr = emg_u_hold(:, tr);
            emg_r_tr = emg_r_hold(:, tr);

            for k = 1:K
                w = tapers(:, k);

                X   = fft(eeg_tr .* w, n_fft);
                Y_u = fft(emg_u_tr .* w, n_fft);
                Y_r = fft(emg_r_tr .* w, n_fft);

                X   = X(1:n_freq_bins);
                Y_u = Y_u(1:n_freq_bins);
                Y_r = Y_r(1:n_freq_bins);

                S_ee        = S_ee + abs(X).^2;
                S_mm_unrect = S_mm_unrect + abs(Y_u).^2;
                S_em_unrect = S_em_unrect + (X .* conj(Y_u));
                S_mm_rect   = S_mm_rect + abs(Y_r).^2;
                S_em_rect   = S_em_rect + (X .* conj(Y_r));
            end
        end

        half_bw = cfg.nw / (n_hold_pnts / fs);
        method_desc = sprintf('DPSS Multitaper (NW = %.1f, K = %d tapers, \\Delta f = \\pm%.2f Hz)', ...
            cfg.nw, K, half_bw);

    case 'welch'
        subwin_pts  = round(cfg.win_sec * fs);
        n_overlap   = round(subwin_pts * (cfg.overlap_pct / 100));
        n_fft       = subwin_pts;
        w_hanning   = hanning(subwin_pts);

        segs_per_trial = floor((n_hold_pnts - n_overlap) / (subwin_pts - n_overlap));
        total_L        = n_valid_trials * segs_per_trial;

        n_freq_bins = (n_fft / 2) + 1;
        S_ee        = zeros(n_freq_bins, 1);
        S_mm_unrect = zeros(n_freq_bins, 1);
        S_em_unrect = zeros(n_freq_bins, 1);
        S_mm_rect   = zeros(n_freq_bins, 1);
        S_em_rect   = zeros(n_freq_bins, 1);

        for tr = 1:n_valid_trials
            [P_ee, f_axis] = cpsd(eeg_hold(:, tr),   eeg_hold(:, tr),   w_hanning, n_overlap, n_fft, fs);
            [P_mm_u, ~]    = cpsd(emg_u_hold(:, tr), emg_u_hold(:, tr), w_hanning, n_overlap, n_fft, fs);
            [P_em_u, ~]    = cpsd(eeg_hold(:, tr),   emg_u_hold(:, tr), w_hanning, n_overlap, n_fft, fs);
            [P_mm_r, ~]    = cpsd(emg_r_hold(:, tr), emg_r_hold(:, tr), w_hanning, n_overlap, n_fft, fs);
            [P_em_r, ~]    = cpsd(eeg_hold(:, tr),   emg_r_hold(:, tr), w_hanning, n_overlap, n_fft, fs);

            S_ee        = S_ee + P_ee;
            S_mm_unrect = S_mm_unrect + P_mm_u;
            S_em_unrect = S_em_unrect + P_em_u;
            S_mm_rect   = S_mm_rect + P_mm_r;
            S_em_rect   = S_em_rect + P_em_r;
        end

        method_desc = sprintf('Welch Periodogram (%.1fs window, %d%% overlap)', ...
            cfg.win_sec, cfg.overlap_pct);
end

cmc_unrect = (abs(S_em_unrect).^2) ./ (S_ee .* S_mm_unrect);
cmc_rect   = (abs(S_em_rect).^2)   ./ (S_ee .* S_mm_rect);

% Theoretical 95% Confidence Limit for 1D spectrum
cl95 = 1 - (0.05)^(1 / (total_L - 1));

% -------------------------------------------------------------------------
% 9. Time-Frequency CMC Estimation (Smoothed Morlet Wavelets)
% -------------------------------------------------------------------------
tf_results = struct();
if cfg.do_tf
    % Downsample time grid for TF estimation (15 ms step)
    time_step_sec = 0.015;
    time_ds       = cfg.epoch_win_sec(1) : time_step_sec : cfg.epoch_win_sec(2);
    n_time_ds     = length(time_ds);
    time_ds_idx   = round((time_ds - cfg.epoch_win_sec(1)) * fs) + 1;
    time_ds_idx   = max(1, min(epoch_pnts, time_ds_idx));

    tf_freqs    = cfg.freq_range(1):1:cfg.freq_range(2);
    n_tf_freqs  = length(tf_freqs);
    n_pnts      = epoch_pnts;
    t_wavelet   = -2:(1/fs):2;
    n_conv      = n_pnts + length(t_wavelet) - 1;
    n_conv_pow2 = 2^nextpow2(n_conv);
    half_w      = (length(t_wavelet) - 1) / 2;

    TF_S_ee        = zeros(n_tf_freqs, n_time_ds);
    TF_S_mm_rect   = zeros(n_tf_freqs, n_time_ds);
    TF_S_em_rect   = zeros(n_tf_freqs, n_time_ds);
    TF_S_mm_unrect = zeros(n_tf_freqs, n_time_ds);
    TF_S_em_unrect = zeros(n_tf_freqs, n_time_ds);

    for fi = 1:n_tf_freqs
        f0 = tf_freqs(fi);
        s0 = cfg.tf_n_cycles / (2 * pi * f0);

        wavelet = exp(2 * 1i * pi * f0 * t_wavelet) .* exp(-t_wavelet.^2 / (2 * s0^2));
        wavelet = wavelet / sum(abs(wavelet));
        fft_w   = fft(wavelet, n_conv_pow2);

        for tr = 1:n_valid_trials
            fft_eeg = fft(eeg_sig(:, tr)', n_conv_pow2);
            fft_mr  = fft(emg_rect(:, tr)', n_conv_pow2);
            fft_mu  = fft(emg_unrect(:, tr)', n_conv_pow2);

            conv_eeg = ifft(fft_eeg .* fft_w, n_conv_pow2);
            conv_mr  = ifft(fft_mr  .* fft_w, n_conv_pow2);
            conv_mu  = ifft(fft_mu  .* fft_w, n_conv_pow2);

            X_full  = conv_eeg(half_w + 1 : half_w + n_pnts);
            Yr_full = conv_mr(half_w + 1 : half_w + n_pnts);
            Yu_full = conv_mu(half_w + 1 : half_w + n_pnts);

            X  = X_full(time_ds_idx);
            Yr = Yr_full(time_ds_idx);
            Yu = Yu_full(time_ds_idx);

            TF_S_ee(fi, :)        = TF_S_ee(fi, :) + (abs(X).^2);
            TF_S_mm_rect(fi, :)   = TF_S_mm_rect(fi, :) + (abs(Yr).^2);
            TF_S_em_rect(fi, :)   = TF_S_em_rect(fi, :) + (X .* conj(Yr));
            TF_S_mm_unrect(fi, :) = TF_S_mm_unrect(fi, :) + (abs(Yu).^2);
            TF_S_em_unrect(fi, :) = TF_S_em_unrect(fi, :) + (X .* conj(Yu));
        end
    end

    % 2D Gaussian Kernel Smoothing across Time and Frequency
    sigma_t_pts = max(1, round(0.08 / time_step_sec)); % ~160 ms window
    sigma_f_pts = 1.2;                                 % ~2.4 Hz window

    % Smooth real and imaginary components independently
    smooth_2d = @(M) complex(imgaussfilt(real(M), [sigma_f_pts, sigma_t_pts]), ...
        imgaussfilt(imag(M), [sigma_f_pts, sigma_t_pts]));

    TF_S_ee_s        = real(smooth_2d(TF_S_ee));
    TF_S_mm_rect_s   = real(smooth_2d(TF_S_mm_rect));
    TF_S_em_rect_s   = smooth_2d(TF_S_em_rect);
    TF_S_mm_unrect_s = real(smooth_2d(TF_S_mm_unrect));
    TF_S_em_unrect_s = smooth_2d(TF_S_em_unrect);

    % Compute Smoothed Coherence
    tf_cmc_rect   = (abs(TF_S_em_rect_s).^2)   ./ (TF_S_ee_s .* TF_S_mm_rect_s);
    tf_cmc_unrect = (abs(TF_S_em_unrect_s).^2) ./ (TF_S_ee_s .* TF_S_mm_unrect_s);

    % Confidence Limits for TF map
    tf_cl95 = 1 - (0.05)^(1 / (n_valid_trials - 1));
    tf_cl99 = 1 - (0.01)^(1 / (n_valid_trials - 1));

    tf_results.time       = time_ds;
    tf_results.freqs      = tf_freqs;
    tf_results.cmc_rect   = tf_cmc_rect;
    tf_results.cmc_unrect = tf_cmc_unrect;
    tf_results.cl95       = tf_cl95;
    tf_results.cl99       = tf_cl99;
end

% -------------------------------------------------------------------------
% 10. Pack Results Struct
% -------------------------------------------------------------------------
cmc_results.freqs          = f_axis;
cmc_results.time_vec       = time_vec;
cmc_results.cmc_unrect     = cmc_unrect;
cmc_results.cmc_rect       = cmc_rect;
cmc_results.cl95           = cl95;
cmc_results.total_L        = total_L;
cmc_results.n_trials       = n_valid_trials;
cmc_results.method         = cfg.method;
cmc_results.eeg_descriptor = eeg_desc;
cmc_results.emg_label      = emg_label;
cmc_results.tf             = tf_results;

% -------------------------------------------------------------------------
% 11. Visualisation (1D Spectrum + Time-Frequency Maps)
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot, return; end

if cfg.do_tf
    fh = figure('Color', 'w', 'Position', [100, 80, 1150, 780], 'Name', 'Corticomuscular Coherence Dashboard');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Panel 1: 1D CMC Spectrum (Top Row, Spanning 2 Columns) ---
    nexttile([1 2]);
    freq_mask   = (f_axis >= cfg.freq_range(1)) & (f_axis <= cfg.freq_range(2));
    f_plot      = f_axis(freq_mask);
    unrect_plot = cmc_unrect(freq_mask);
    rect_plot   = cmc_rect(freq_mask);

    y_max = max([max(unrect_plot), max(rect_plot), cl95 * 1.5, 0.05]) * 1.15;

    patch([15 30 30 15], [0 0 y_max y_max], [0.92 0.95 0.98], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'Beta Band (15-30 Hz)');
    hold on;
    patch([30 45 45 30], [0 0 y_max y_max], [0.98 0.94 0.92], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Gamma Band (30-45 Hz)');

    h1 = plot(f_plot, rect_plot, 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, ...
        'DisplayName', 'Rectified EMG (|EMG|)');
    h2 = plot(f_plot, unrect_plot, 'Color', [0.15 0.45 0.85], 'LineWidth', 1.8, ...
        'DisplayName', 'Unrectified EMG (Raw)');
    h_cl = plot([cfg.freq_range(1), cfg.freq_range(2)], [cl95, cl95], '--k', ...
        'LineWidth', 1.3, 'DisplayName', sprintf('95%% Conf. Limit (CL_{95} = %.4f, L = %d)', cl95, total_L));

    xlim(cfg.freq_range);
    ylim([0, y_max]);
    xlabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Coherence', 'FontSize', 10, 'FontWeight', 'bold');
    title({sprintf('Stationary Task CMC (+%.1fs to +%.1fs post [%s]): %s vs. %s', ...
        cfg.task_win_sec(1), cfg.task_win_sec(2), string(cfg.event_marker), eeg_desc, emg_label), method_desc}, ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    legend([h1, h2, h_cl], 'Location', 'northeast', 'FontSize', 8);

    % Common colour scaling for TF plots
    max_tf_val = max([max(tf_cmc_rect(:)), max(tf_cmc_unrect(:)), tf_cl95 * 1.5, 0.05]);

    % --- Panel 2: TF Map for Rectified EMG (Bottom-Left) ---
    nexttile;
    contourf(time_ds, tf_freqs, tf_cmc_rect, 40, 'LineColor', 'none');
    hold on;
    contour(time_ds, tf_freqs, tf_cmc_rect, [tf_cl95, tf_cl95], 'LineColor', 'k', 'LineWidth', 1.4);
    xline(0.0, '--w', 'Cue (21)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    xline(5.0, '--w', 'Release (22)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    colormap(gca, parula(256));
    caxis([0, max_tf_val]);
    colorbar;
    xlim(cfg.epoch_win_sec);
    xlabel('Time relative to stimulus ''21'' (s)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('TF Coherence: Rectified |EMG| (Black contour = CL_{95} [%.3f])', tf_cl95), ...
        'FontSize', 10, 'FontWeight', 'bold');
    grid on;

    % --- Panel 3: TF Map for Unrectified EMG (Bottom-Right) ---
    nexttile;
    contourf(time_ds, tf_freqs, tf_cmc_unrect, 40, 'LineColor', 'none');
    hold on;
    contour(time_ds, tf_freqs, tf_cmc_unrect, [tf_cl95, tf_cl95], 'LineColor', 'k', 'LineWidth', 1.4);
    xline(0.0, '--w', 'Cue (21)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    xline(5.0, '--w', 'Release (22)', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    colormap(gca, parula(256));
    caxis([0, max_tf_val]);
    colorbar;
    xlim(cfg.epoch_win_sec);
    xlabel('Time relative to stimulus ''21'' (s)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('TF Coherence: Unrectified Raw EMG (Black contour = CL_{95} [%.3f])', tf_cl95), ...
        'FontSize', 10, 'FontWeight', 'bold');
    grid on;

else
    % Standard 1D plot only
    fh = figure('Color', 'w', 'Position', [150, 150, 950, 550], 'Name', 'Task-Locked CMC');
    freq_mask   = (f_axis >= cfg.freq_range(1)) & (f_axis <= cfg.freq_range(2));
    f_plot      = f_axis(freq_mask);
    unrect_plot = cmc_unrect(freq_mask);
    rect_plot   = cmc_rect(freq_mask);

    y_max = max([max(unrect_plot), max(rect_plot), cl95 * 1.5, 0.05]) * 1.15;

    patch([15 30 30 15], [0 0 y_max y_max], [0.92 0.95 0.98], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    hold on;
    patch([30 45 45 30], [0 0 y_max y_max], [0.98 0.94 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    h1 = plot(f_plot, rect_plot, 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'DisplayName', 'Rectified EMG (|EMG|)');
    h2 = plot(f_plot, unrect_plot, 'Color', [0.15 0.45 0.85], 'LineWidth', 1.8, 'DisplayName', 'Unrectified EMG (Raw)');
    h_cl = plot([cfg.freq_range(1), cfg.freq_range(2)], [cl95, cl95], '--k', 'LineWidth', 1.3, ...
        'DisplayName', sprintf('95%% Conf. Limit (CL_{95} = %.4f)', cl95));

    xlim(cfg.freq_range);
    ylim([0, y_max]);
    xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Magnitude-Squared Coherence', 'FontSize', 11, 'FontWeight', 'bold');
    title({sprintf('Stationary Task CMC (+%.1fs to +%.1fs post [%s]): %s vs. %s', ...
        cfg.task_win_sec(1), cfg.task_win_sec(2), string(cfg.event_marker), eeg_desc, emg_label), method_desc}, ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    legend([h1, h2, h_cl], 'Location', 'northeast', 'FontSize', 9);
end

% -------------------------------------------------------------------------
% 12. Console Summary
% -------------------------------------------------------------------------
beta_mask = (f_plot >= 15) & (f_plot <= 30);
[max_beta_rect, idx_b_rect] = max(rect_plot(beta_mask));
[max_beta_unrect, idx_b_un] = max(unrect_plot(beta_mask));
f_beta = f_plot(beta_mask);

fprintf('\n===================================================================\n');
fprintf('  Task-Locked Corticomuscular Coherence Summary\n');
fprintf('===================================================================\n');
fprintf('  Spectral Method:           %s\n', method_desc);
fprintf('  Epoch Window:              [%.1fs, +%.1fs] relative to [%s]\n', ...
    cfg.epoch_win_sec(1), cfg.epoch_win_sec(2), string(cfg.event_marker));
fprintf('  1D Hold Window:            [+%.1fs, +%.1fs]\n', cfg.task_win_sec(1), cfg.task_win_sec(2));
fprintf('  Valid Task Epochs:         N = %d (Total Estimates L = %d)\n', n_valid_trials, total_L);
fprintf('  EEG Derivation:            %s\n', eeg_desc);
fprintf('  EMG Channel:               %s (Highpass: %d Hz)\n', emg_label, cfg.emg_hp);
fprintf('  95%% Confidence Threshold:  %.4f\n', cl95);
fprintf('  -----------------------------------------------------------------\n');
fprintf('  Peak Beta CMC (Rectified):   %.4f at %.2f Hz %s\n', ...
    max_beta_rect, f_beta(idx_b_rect), get_sig_label(max_beta_rect, cl95));
fprintf('  Peak Beta CMC (Unrectified): %.4f at %.2f Hz %s\n', ...
    max_beta_unrect, f_beta(idx_b_un), get_sig_label(max_beta_unrect, cl95));
fprintf('===================================================================\n\n');

end

function lbl = get_sig_label(val, cl)
if val >= cl
    lbl = '[SIGNIFICANT]';
else
    lbl = '[Non-Significant]';
end
end