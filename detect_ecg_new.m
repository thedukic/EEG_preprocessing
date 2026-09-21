function [ECGmask, ECGEpochs, ECGlatency, heartData, brainData, pulsEstimate] = detect_ecg_new(DATA, cfg)
% DETECT_ECG Detects R-peaks and extracts QRS epochs from continuous ECG data.
% Uses local dynamic thresholding and gap search-back for robust detection
% across low-amplitude and respiration-modulated BioSemi recordings.
%
% Inputs:
%    DATA : EEGLAB dataset structure containing the ECG channel
%    cfg  : Configuration structure with optional fields:
%           - cfg.win_heart     : [start_ms, end_ms] relative to R-peak (default: [-200, 250])
%           - cfg.plot_path     : Directory path to save figure (default: DATA.ALSUTRECHT.subject.figures)
%           - cfg.plot_visible  : 'on' or 'off' (default: 'off')
%           - cfg.do_plot       : true/false whether to generate the plot (default: true)
%
% Outputs:
%    ECGmask      : Logical mask across continuous samples indicating QRS periods
%    ECGEpochs    : [N_epochs x 2] Start and end sample indices of valid epochs
%    ECGlatency   : [1 x N_epochs] Time latency (seconds) of detected R-peaks
%    heartData    : Continuous raw ECG channel data (polarity-corrected)
%    brainData    : Empty array for pipeline compatibility
%    pulsEstimate : Estimated median heart rate in beats per minute (bpm)

%% 1. Parameter Parsing and Initialisation
ECGmask      = NaN;
ECGEpochs    = NaN;
ECGlatency   = NaN;
heartData    = NaN;
brainData    = [];
pulsEstimate = NaN;

if nargin < 2, cfg = struct(); end

% Support backwards compatibility if numeric window is passed directly
if isnumeric(cfg)
    tmp_win = cfg;
    cfg = struct();
    cfg.win_heart = tmp_win;
end

% Parse options with robust fallbacks
if isscalar(cfg.win_heart)
    cfg.win_heart = [-abs(cfg.win_heart), abs(cfg.win_heart)];
end

if ~isfield(cfg, 'do_plot') || isempty(cfg.do_plot)
    cfg.do_plot = true;
end

% Validate ECG lead
channel_ecg = strcmp({DATA.chanlocs.labels}, 'ECG');
if ~any(channel_ecg) || sum(channel_ecg) ~= 1
    warning('ECG channel not found or multiple matches present.');
    return;
end

fs = DATA.srate;
heartData = double(DATA.data(channel_ecg, :));
n_samples = length(heartData);

%% 2. Dedicated QRS Filtering (6 to 18 Hz)
% % Suppresses slow baseline drift (<6 Hz), cortical alpha, and high-frequency EMG (>18 Hz)
% [b_bp, a_bp] = butter(3, [6, 18] / (fs / 2), 'bandpass');
% ecg_filt = filtfilt(b_bp, a_bp, heartData);
is_approximated = isfield(DATA, 'ALSUTRECHT') && ...
    isfield(DATA.ALSUTRECHT, 'subject') && ...
    isfield(DATA.ALSUTRECHT.subject, 'ecg') && ...
    contains(lower(DATA.ALSUTRECHT.subject.ecg), 'approx');

if is_approximated
    % Shift above occipital alpha (10 Hz) and below muscle beta
    [b_bp, a_bp] = butter(3, [14, 28] / (fs / 2), 'bandpass');
else
    % Standard thoracic lead: wider cardiac power band
    [b_bp, a_bp] = butter(3, [6, 18] / (fs / 2), 'bandpass');
end
ecg_filt = filtfilt(b_bp, a_bp, heartData);

%% 3. Energy Envelope (Derivative + Moving Average Integration)
ecg_diff = diff([ecg_filt(1), ecg_filt]);
ecg_sq   = ecg_diff.^2;

int_samples = round(0.090 * fs); % 90 ms integration window
b_ma = (1 / int_samples) * ones(1, int_samples);
y_energy = filtfilt(b_ma, 1, ecg_sq);

% Clamp extreme motion outliers to prevent local gain suppression
med_eng = median(y_energy);
p95_eng = prctile(y_energy, 95);
cap_val = med_eng + 15 * max(p95_eng - med_eng, eps);
y_energy_capped = min(y_energy, cap_val);

%% 4. Local Dynamic Thresholding (Automatic Gain Tracking)
% Continuous sliding 4.0-second window tracks local baseline and peak amplitude
win_track  = round(4.0 * fs);
local_base = movmedian(y_energy_capped, win_track);
local_max  = movmax(y_energy_capped, win_track);

% Dynamic threshold: baseline + 20% of local peak excursion (with floor protection)
min_excursion = 0.15 * max(med_eng, eps);
local_thresh  = local_base + max(0.20 * (local_max - local_base), min_excursion);

% Normalise energy against the dynamic threshold (true peaks cross 1.0)
y_norm = y_energy_capped ./ local_thresh;

%% 5. Candidate Peak Detection with Search-Back
min_dist_samples = round(0.350 * fs); % Max physiological rate ~171 bpm

[~, peak_locs] = findpeaks(y_norm, ...
    'MinPeakHeight', 1.0, ...
    'MinPeakDistance', min_dist_samples);

% Pan-Tompkins search-back: recover missed beats in unusually wide RR gaps
if length(peak_locs) >= 4
    rr_init = diff(peak_locs);
    valid_rr = rr_init(rr_init >= round(0.40 * fs) & rr_init <= round(1.80 * fs));

    if ~isempty(valid_rr)
        med_rr = median(valid_rr);
        recovered_locs = [];

        for i_gap = 1:length(peak_locs) - 1
            gap_len = peak_locs(i_gap + 1) - peak_locs(i_gap);
            if gap_len > 1.60 * med_rr
                % Search within the gap using a relaxed threshold (0.45)
                g_start = peak_locs(i_gap) + round(0.25 * fs);
                g_end   = peak_locs(i_gap + 1) - round(0.25 * fs);

                if g_end > g_start
                    [~, gap_pks] = findpeaks(y_norm(g_start:g_end), ...
                        'MinPeakHeight', 0.45, ...
                        'MinPeakDistance', min_dist_samples);
                    if ~isempty(gap_pks)
                        recovered_locs = [recovered_locs, g_start + gap_pks(:)' - 1]; %#ok<AGROW>
                    end
                end
            end
        end

        if ~isempty(recovered_locs)
            peak_locs = sort([peak_locs(:); recovered_locs(:)]);
        end
    end
end

if isempty(peak_locs) || length(peak_locs) < 5
    warning('ECG signal contains insufficient detectable peaks. Quitting...');
    return;
end

%% 6. Local Polarity Consensus
test_win = round(0.040 * fs);
peak_amps = zeros(length(peak_locs), 1);
for k = 1:length(peak_locs)
    idx1 = max(1, peak_locs(k) - test_win);
    idx2 = min(n_samples, peak_locs(k) + test_win);
    [~, max_rel] = max(abs(ecg_filt(idx1:idx2)));
    peak_amps(k) = ecg_filt(idx1 + max_rel - 1);
end

if median(peak_amps) < 0
    fprintf('Inverted ECG detected via peak consensus. Flipping signal.\n');
    heartData = -heartData;
    ecg_filt  = -ecg_filt;
end

%% 7. Exact R-Peak Alignment
search_win = round(0.060 * fs);
true_locs = zeros(size(peak_locs));

for k = 1:length(peak_locs)
    idx_start = max(1, peak_locs(k) - search_win);
    idx_end   = min(n_samples, peak_locs(k) + search_win);

    [~, local_max] = max(heartData(idx_start:idx_end));
    true_locs(k) = idx_start + local_max - 1;
end
true_locs = unique(true_locs);
n_raw_peaks = length(true_locs);

%% 8. Epoch Extraction and Boundary Checking
win_samples = round(abs(cfg.win_heart) ./ (1000 / fs));

starts = true_locs - win_samples(1);
ends   = true_locs + win_samples(2);

valid_bounds = (starts >= 1) & (ends <= n_samples);
starts    = starts(valid_bounds);
ends      = ends(valid_bounds);
true_locs = true_locs(valid_bounds);

n_epochs  = length(starts);
epoch_len = win_samples(1) + win_samples(2) + 1;
ECG = zeros(n_epochs, epoch_len);

for k = 1:n_epochs
    ECG(k, :) = heartData(starts(k):ends(k));
end

ECGEpochs  = [starts(:), ends(:)];
ECGlatency = DATA.times(true_locs) / 1000;

%% 9. Outlier Rejection via Template Matching
T = linspace(cfg.win_heart(1), cfg.win_heart(2), epoch_len) / 1000;

[valid_idx, ~, ~] = find_template_outliers(ECG', T);
ECGEpochs  = ECGEpochs(valid_idx, :);
ECG        = ECG(valid_idx, :);
ECGlatency = ECGlatency(valid_idx);
n_epochs   = size(ECG, 1);

fprintf('Peak detection summary: %d raw peaks detected -> %d retained after outlier rejection (%d rejected).\n', ...
    n_raw_peaks, n_epochs, n_raw_peaks - n_epochs);

if n_epochs == 0
    warning('All detected ECG epochs were rejected as outliers.');
    return;
end

% Build binary continuous mask
ECGmask = false(1, n_samples);
for k = 1:n_epochs
    ECGmask(ECGEpochs(k, 1):ECGEpochs(k, 2)) = true;
end

%% 10. Heart Rate Estimation
ibi = diff(ECGlatency);
valid_ibi = ibi(ibi >= 0.40 & ibi <= 2.0); % 30 to 150 bpm
if ~isempty(valid_ibi)
    pulsEstimate = 60 / median(valid_ibi);
end

fprintf('Final valid QRS complexes: %d (Estimated rate: %.1f bpm).\n', n_epochs, pulsEstimate);

%% 11. Diagnostics Plotting
if ~cfg.do_plot
    return;
end

fh = figure('Color', 'w', 'Position', [100, 100, 800, 450], 'Visible', cfg.plot_visible);
hold on;

colors = brewermap(n_epochs, 'YlGnBu');
for k = 1:n_epochs
    h_line = plot(T, ECG(k, :), 'LineWidth', 0.8, 'Color', colors(k, :));
    h_line.Color(4) = 0.12;
end

plot(T, mean(ECG, 1), 'Color', [0.65, 0.10, 0.10], 'LineWidth', 2.5);
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Box', 'off', 'FontName', 'Helvetica');
xlabel('Time Relative to R-Peak (s)', 'FontWeight', 'bold');
ylabel('Amplitude (\muV)', 'FontWeight', 'bold');
title(sprintf('Detected QRS (N = %d, Rate = %.1f bpm)\nfrom %s ECG', n_epochs, pulsEstimate, DATA.ALSUTRECHT.subject.ecg), 'FontWeight', 'bold');
hold off;

% Save figure
sub_id = DATA.ALSUTRECHT.subject.id;
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [sub_id '_detected_ecg'], [20, 11]);

end