function smoothed_psd = detrend_spectra(psdspectra, freq)

num_channel = size(psdspectra, 1);
smoothed_psd = NaN(size(psdspectra));
log_freq = log10(freq(:));

for i_channel = 1:num_channel
    % --- AUTOMATIC PEAK DETECTION LOGIC ---
    % 1. Calculate mean spectrum and log-transform
    log_psd  = log10(psdspectra(i_channel, :))';

    % 2. Simple 1/f fitting (linear fit in log-log space)
    fit_idx = freq > 1 & freq < 60;
    cen_idx = freq > 4 & freq < 40;
    fit_idx(cen_idx) = false;

    p = polyfit(log_freq(fit_idx), log_psd(fit_idx), 1);
    aperiodic_component = p(1) * log_freq + p(2);

    % 3. Flatten the spectrum (residual)
    flattened_psd = log_psd - aperiodic_component;

    % Smooth slightly to remove the high-frequency "wobbles" (false positives)
    smoothed_psd(i_channel, :) = smoothdata(flattened_psd, 'gaussian', 5);

    % % =========================================================
    % % --- NEW: REMOVE BROAD EMG "BELLY" (ROLLING BASELINE) ---
    % % =========================================================
    % % Calculate frequency resolution to define window size in bins
    % df = freq(2) - freq(1);
    % 
    % % Define a window wide enough to ignore sharp peaks,
    % % but narrow enough to track the broad belly (e.g., 4 Hz)
    % window_hz = 4;
    % window_bins = round(window_hz / df);
    % 
    % % Calculate the rolling baseline using a moving median
    % broad_baseline(i_channel, :) = movmedian(smoothed_psd, window_bins);

    % % Subtract the belly from the smoothed spectrum
    % detrended_psd(i_channel, :) = smoothed_psd - broad_baseline;
end

end