function [potential_freqs, fh] = find_peaks(psdspectra, freq, optVisible)

% Define
do_plot = true;
idx_freq = freq > 0.5;
psdspectra = psdspectra(idx_freq,:);
freq = freq(idx_freq);

% --- AUTOMATIC PEAK DETECTION LOGIC ---
% 1. Calculate mean spectrum and log-transform
% Recommendation: If peaks are highly focal (1-2 bad electrodes),
% prctile(psdspectra, 95, 2) is often much better than mean()
mean_psd = mean(psdspectra, 2);

log_freq = log10(freq(:));
log_psd  = log10(mean_psd(:));

% 2. Simple 1/f fitting (linear fit in log-log space)
fit_idx = freq > 0.5 & freq < 70;
cen_idx = freq > 1.5 & freq < 35;
fit_idx(cen_idx) = false;

p = polyfit(log_freq(fit_idx), log_psd(fit_idx), 1);
aperiodic_component = p(1) * log_freq + p(2);

% 3. Flatten the spectrum (residual)
flattened_psd = log_psd - aperiodic_component;

% Smooth slightly to remove the high-frequency "wobbles" (false positives)
smoothed_psd = smoothdata(flattened_psd, 'gaussian', 5);

% =========================================================
% --- NEW: REMOVE BROAD EMG "BELLY" (ROLLING BASELINE) ---
% =========================================================
% Calculate frequency resolution to define window size in bins
df = freq(2) - freq(1);

% Define a window wide enough to ignore sharp peaks,
% but narrow enough to track the broad belly (e.g., 4 Hz)
window_hz = 4;
window_bins = round(window_hz / df);

% Calculate the rolling baseline using a moving median
broad_baseline = movmedian(smoothed_psd, window_bins);

% Subtract the belly from the smoothed spectrum
detrended_psd = smoothed_psd - broad_baseline;
% =========================================================

% --- DYNAMIC NOISE THRESHOLDING ---
% Now calculate the noise floor on the DETRENDED spectrum
hf_idx = freq > 40;
hf_data = detrended_psd(hf_idx);

% Median Absolute Deviation (MAD)
noise_mad = median(abs(hf_data - median(hf_data)));
noise_std = noise_mad * 1.4826;

sigma_multiplier = 5;
dynamic_prom = max(0.1, noise_std * sigma_multiplier);

fprintf('Baseline High-Freq Noise (Std): %.3f | Dynamic Prominence: %.3f\n', noise_std, dynamic_prom);

% 4. Detect narrow peaks (Using the detrended PSD)
[pks, candidate_locs, widths, proms] = findpeaks(detrended_psd, freq, ...
    'MinPeakProminence', dynamic_prom, ...
    'MaxPeakWidth', 1);

% 5. Identify high-frequency peaks (40 to max Hz)
fundamental_range = [40, max(freq)];
sel_idx1 = candidate_locs >= fundamental_range(1) & candidate_locs <= fundamental_range(2);
sel_idx2 = pks > dynamic_prom;
sel_idx = sel_idx1 & sel_idx2;

potential_freqs = candidate_locs(sel_idx);
potential_pks = pks(sel_idx);

if do_plot % && ~isempty(potential_freqs)
    fh = figure('Color', 'w', 'Name', 'Peak Detection Debugger', 'Visible', optVisible);

    % Subplot 1: The Fit
    subplot(2,1,1); hold on;
    plot(freq, log_psd, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original (Log)');
    plot(freq, aperiodic_component, 'r--', 'LineWidth', 1.5, 'DisplayName', '1/f Fit');
    grid on; xlabel('Frequency (Hz)'); ylabel('log_{10} Power');
    title('1/f Background Fit');
    legend('Location', 'northeast');
    xlim([0 100]);

    % Subplot 2: The Flattened Spectrum & Found Peaks
    subplot(2,1,2); hold on;

    % Plot the raw flattened spectrum (with the belly) in light blue
    plot(freq, smoothed_psd, 'Color', [0.6 0.8 1], 'DisplayName', 'Smoothed Spectrum (w/ Belly)');

    % Plot the moving median baseline (the tracked belly)
    plot(freq, broad_baseline, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Rolling Baseline (EMG Belly)');

    % Plot the final detrended spectrum that findpeaks uses
    plot(freq, detrended_psd, 'b', 'LineWidth', 1.5, 'DisplayName', 'Detrended (Flattened) Spectrum');

    % Plot all candidate peaks
    scatter(candidate_locs, pks, 50, [0.3 0.3 0.3], 'DisplayName', 'All Narrow Peaks');

    % Highlight the identified target peaks in red
    scatter(potential_freqs, potential_pks, 50, 'r', 'filled', 'DisplayName', 'Identified High-Freq Peaks');

    % Add a visual representation of the dynamic threshold
    yline(dynamic_prom, 'k--', sprintf('Prominence Threshold (%.2f)', dynamic_prom), 'LabelHorizontalAlignment', 'left');

    grid on; xlabel('Frequency (Hz)'); ylabel('Relative Power (dB-like)');
    title(sprintf('Peak Detection (Dynamic Prominence: %.3f)', dynamic_prom));
    legend('Location', 'northeast');
    xlim([0 100]);
end

if length(potential_freqs) < 3
    fprintf('Too few peaks found in high-frequency range (N = %d).\n', length(potential_freqs));
    potential_freqs = [];
end
end