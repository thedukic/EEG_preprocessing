% =========================================================================
% Whole-Head CMC: 128 Channels (Unrectified vs Rectified EMG)
% =========================================================================
[EEG, EMG, EXT] = separate_electrodetypes(DATA(1));





% 1. Define parameters and channels
n_eeg_ch     = 128;
eeg_channels = 1:n_eeg_ch;
emg_ch_idx   = [129 130];
emg_ch_idx   = [131 132];
fs           = 2048;
target_event = 'condition 21';

% Define Beta band for topoplot
beta_range = [15 30]; % Hz

% 2. High-pass filter continuous data (Zero-phase)
% EEG: 1 Hz high-pass to remove slow drifts
[b_eeg, a_eeg] = butter(2, 1 / (fs / 2), 'high');
% EMG: 10 Hz high-pass (crucial before rectification)
[b_emg, a_emg] = butter(2, 10 / (fs / 2), 'high');

% Apply to EEG channels
EEG_filtered = zeros(n_eeg_ch, size(EEG.data, 2));
for ch = 1:n_eeg_ch
    EEG_filtered(ch, :) = filtfilt(b_eeg, a_eeg, double(EEG.data(ch, :)));
end

% Extract and filter continuous bipolar EMG
EMG_bipolar_raw = double(EEG.data(emg_ch_idx(1), :)) - double(EEG.data(emg_ch_idx(2), :));
EMG_bipolar     = filtfilt(b_emg, a_emg, EMG_bipolar_raw);

% 3. Common Average Reference (CAR) for filtered EEG
eeg_mean = mean(EEG_filtered, 1);
EEG_data_car = EEG_filtered - eeg_mean;

% 4. Spectral parameters
win_len  = fs; % 1-second window (1 Hz resolution)
window   = hann(win_len);
noverlap = floor(0.5 * win_len);
nfft     = win_len;
n_freqs  = floor(nfft/2) + 1; % Number of frequency bins returned by pwelch

% Pre-allocate spectral sums (freqs x channels)
Pxx_sum = zeros(n_freqs, n_eeg_ch);
Pyy_unrect_sum = zeros(n_freqs, 1);
Pxy_unrect_sum = zeros(n_freqs, n_eeg_ch);
Pyy_rect_sum = zeros(n_freqs, 1);
Pxy_rect_sum = zeros(n_freqs, n_eeg_ch);

% 5. Epoching and Spectral Calculation Loop
all_events = strtrim(string({EEG.event.type}));
target_indices = find(strcmpi(all_events, target_event) | strcmpi(all_events, '21'));
n_trials = length(target_indices);

if n_trials == 0
    error('No events found matching %s.', target_event);
end

% Define window offsets in samples (+1s to +5s)
offset_start = 1 * fs;
offset_end   = 5 * fs - 1;
n_samples = offset_end - offset_start + 1;

valid_trials = 0;

for tr = 1:n_trials
    ev_latency = round(EEG.event(target_indices(tr)).latency);
    start_samp = ev_latency + offset_start;
    end_samp   = ev_latency + offset_end;

    if start_samp < 1 || end_samp > size(EEG_data_car, 2)
        continue;
    end

    % Extract trial matrices: EEG is [samples x channels]
    x = EEG_data_car(:, start_samp:end_samp)';
    y_unrect = EMG_bipolar(start_samp:end_samp)';

    % Rectify EMG (now safely high-pass filtered)
    y_rect = abs(y_unrect);

    % Remove the mean (DC offset) for the local epoch
    x = x - mean(x, 1);
    y_unrect = y_unrect - mean(y_unrect);
    y_rect = y_rect - mean(y_rect);

    % Compute auto-spectra
    [pxx, f] = pwelch(x, window, noverlap, nfft, fs);
    [pyy_u, ~] = pwelch(y_unrect, window, noverlap, nfft, fs);
    [pyy_r, ~] = pwelch(y_rect, window, noverlap, nfft, fs);

    % Replicate EMG vector to match EEG matrix size for vectorized cpsd
    Y_u = repmat(y_unrect, 1, n_eeg_ch);
    Y_r = repmat(y_rect, 1, n_eeg_ch);

    % Compute cross-spectra
    [pxy_u, ~] = cpsd(x, Y_u, window, noverlap, nfft, fs);
    [pxy_r, ~] = cpsd(x, Y_r, window, noverlap, nfft, fs);

    % Accumulate sums
    Pxx_sum = Pxx_sum + pxx;
    Pyy_unrect_sum = Pyy_unrect_sum + pyy_u;
    Pxy_unrect_sum = Pxy_unrect_sum + pxy_u;

    Pyy_rect_sum = Pyy_rect_sum + pyy_r;
    Pxy_rect_sum = Pxy_rect_sum + pxy_r;

    valid_trials = valid_trials + 1;
end

% 6. Compute Magnitude Squared Coherence (MSC)
% Broadcast Pyy (freqs x 1) to match Pxx (freqs x channels)
Pyy_unrect_mat = repmat(Pyy_unrect_sum, 1, n_eeg_ch);
Pyy_rect_mat   = repmat(Pyy_rect_sum, 1, n_eeg_ch);

coh_unrect = (abs(Pxy_unrect_sum).^2) ./ (Pxx_sum .* Pyy_unrect_mat);
coh_rect   = (abs(Pxy_rect_sum).^2)   ./ (Pxx_sum .* Pyy_rect_mat);

% Transpose to output requested [channels x freqs] format
coh_unrect = coh_unrect';
coh_rect   = coh_rect';

% 7. Calculate 95% Confidence Limit
n_segments_per_trial = floor((n_samples - noverlap) / (win_len - noverlap));
L_total = valid_trials * n_segments_per_trial;
alpha = 0.05;
conf_limit = 1 - (alpha)^(1 / (L_total - 1));

% 8. Extract Beta Band Coherence for Topoplots
beta_idx = find(f >= beta_range(1) & f <= beta_range(2));
beta_coh_unrect = mean(coh_unrect(:, beta_idx), 2);
beta_coh_rect   = mean(coh_rect(:, beta_idx), 2);

% 9. Visualisation
figure('Color', 'w', 'Position', [100, 100, 1200, 800], ...
    'Name', 'Whole-Head CMC: Unrectified vs Rectified EMG');

freq_limit_idx = find(f <= 50, 1, 'last'); % Plot up to 50 Hz

% Set a common color scale for both heatmaps based on the max value found
cmax = max([max(coh_unrect(:, 1:freq_limit_idx), [], 'all'), ...
    max(coh_rect(:, 1:freq_limit_idx), [], 'all')]) * 1.1;

% --- Panel 1: Heatmap (Unrectified) ---
subplot(2, 2, 1);
imagesc(f(1:freq_limit_idx), 1:n_eeg_ch, coh_unrect(:, 1:freq_limit_idx));
set(gca, 'YDir', 'normal');
colormap(gca, 'jet');
clim([0 cmax]);
colorbar;
xlabel('Frequency (Hz)');
ylabel('EEG Channel');
title('Unrectified EMG Coherence (1-128)');

% --- Panel 2: Heatmap (Rectified) ---
subplot(2, 2, 2);
imagesc(f(1:freq_limit_idx), 1:n_eeg_ch, coh_rect(:, 1:freq_limit_idx));
set(gca, 'YDir', 'normal');
colormap(gca, 'jet');
clim([0 cmax]);
colorbar;
xlabel('Frequency (Hz)');
ylabel('EEG Channel');
title('Rectified EMG Coherence (1-128)');

% Set a common color scale for topoplots
cmax_topo = max([max(beta_coh_unrect), max(beta_coh_rect)]);

% --- Panel 3: Topoplot (Unrectified) ---
sph = subplot(2, 2, 3);
mytopoplot(beta_coh_unrect, [], sprintf('Unrectified Beta CMC (%d-%d Hz)', beta_range(1), beta_range(2)), sph, [0 cmax_topo]);
colorbar;

% --- Panel 4: Topoplot (Rectified) ---
sph = subplot(2, 2, 4);
mytopoplot(beta_coh_rect, [], sprintf('Rectified Beta CMC (%d-%d Hz)', beta_range(1), beta_range(2)), sph, [0 cmax_topo]);
colorbar;