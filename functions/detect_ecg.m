function [ECGmask, ECGEpochs, ECGlatency, heartData, brainData, pulsEstimate] = detect_ecg(DATA, winHeart, optVisible)

% Initialise
ECGmask       = NaN;
ECGEpochs     = NaN;
ECGlatency    = NaN;
heartData     = NaN;
brainData     = NaN;
pulsEstimate  = NaN;

% Select ECG
ECGtype = DATA.ALSUTRECHT.subject.ecg;
assert(ischar(ECGtype) & ~strcmpi(ECGtype, 'none'));
channel_ecg = strcmp({DATA.chanlocs.labels}, 'ECG');

% Was ECG recorded?
if any(channel_ecg) && sum(channel_ecg) == 1
    heartData = DATA.data(channel_ecg, :);
    brainData = [];

    % MATLAB:
    % https://nl.mathworks.com/help/wavelet/ug/r-wave-detection-in-the-ecg.html
    wt = modwt(heartData, 5);
    wtrec = zeros(size(wt));
    wtrec(4:5, :) = wt(4:5, :);
    y = imodwt(wtrec, 'sym4');

    % Detect and flip inverted ECG signals using the filtered signal 'y'
    % This prevents raw baseline drifts from causing false inversions
    if skewness(y) < 0
        fprintf('Inverted ECG signal detected. Flipping the signal.\n');
        heartData = -heartData;
        y = -y;
    end

    % ---------------------------------------------------------------------
    % FIX 1: Strict QRS Bandpass Filtering
    % Aggressively attenuate T-waves (< 10 Hz) and high-frequency noise (> 20 Hz)
    [b_bp, a_bp] = butter(4, [10 20]/(DATA.srate/2), 'bandpass');
    y_filtered = filtfilt(b_bp, a_bp, y);

    % Accentuate the ECG peaks using the strictly filtered signal
    y_diff = diff([y_filtered(1) y_filtered]);
    y_sq = y_diff.^2;

    % ---------------------------------------------------------------------
    % FIX 2: Narrower Integration Window
    % Shrink the window to 80 ms to strictly capture the narrow QRS energy
    % without accumulating the broader T-wave energy.
    windowLength = round(0.080 * DATA.srate);
    b = (1/windowLength) * ones(1, windowLength);
    a = 1;
    y_energy = filtfilt(b, a, y_sq);

    % ---------------------------------------------------------------------
    % Robust Thresholding for Energy Signals
    [tmp_peaks, ~] = findpeaks(y_energy, 'MinPeakHeight', mean(y_energy), 'MinPeakDistance', 0.25 * DATA.srate);

    if ~isempty(tmp_peaks)
        threshold = 0.5 * median(tmp_peaks);
    else
        threshold = mean(y_energy);
    end

    if isscalar(winHeart)
        fprintf('Heartbeats (L = +-%d ms) detected using a 50%% median peak threshold = %1.2f\n', winHeart, threshold);
    else
        fprintf('Heartbeats (L = %d - %d ms) detected using a 50%% median peak threshold = %1.2f\n', winHeart, threshold);
    end

    % ---------------------------------------------------------------------
    % Detect Peaks (Two-Stage Method)
    [~, approx_locs_samples] = findpeaks(y_energy, 'MinPeakHeight', threshold, 'MinPeakDistance', 0.25 * DATA.srate);

    if isempty(approx_locs_samples)
        warning('Low quality ECG data. Quitting...');
        return;
    end

    search_window = round(0.050 * DATA.srate);
    true_locs_samples = zeros(size(approx_locs_samples));
    qrspeaks = zeros(size(approx_locs_samples));

    for i = 1:length(approx_locs_samples)
        idx_start = max(1, approx_locs_samples(i) - search_window);
        idx_end   = min(length(heartData), approx_locs_samples(i) + search_window);

        [peak_val, local_idx] = max(heartData(idx_start:idx_end));
        true_locs_samples(i) = idx_start + local_idx - 1;
        qrspeaks(i) = peak_val;
    end

    times = DATA.times / 1000;
    ECGlatency = times(true_locs_samples);

    fprintf('Detected %d QRS peaks in %1.1f min of data.\n', length(qrspeaks), times(end)/60);

    % ---------------------------------------------------------------------
    % Select the ECG window
    if isscalar(winHeart)
        winHeart = winHeart * [-1 1];
    end
    assert(length(winHeart) == 2);

    winECGsamples = round(abs(winHeart) ./ (1000/DATA.srate));
    starts = true_locs_samples - winECGsamples(1);
    ends   = true_locs_samples + winECGsamples(2);

    % Identify epochs that fit entirely within the data bounds
    valid_idx = (starts >= 1) & (ends <= DATA.pnts);
    fprintf('Detections removed due to boundary exceedance: %d\n', sum(~valid_idx));

    starts = starts(valid_idx);
    ends   = ends(valid_idx);
    ECGlatency = ECGlatency(valid_idx);

    % Extract valid epochs
    num_time = winECGsamples(1) + winECGsamples(2) + 1;
    num_epoch = length(starts);

    ECG = zeros(num_epoch, num_time);
    for i = 1:num_epoch
        ECG(i, :) = heartData(starts(i):ends(i));
    end

    % ECG epochs
    ECGEpochs = [starts', ends'];

    % ---------------------------------------------------------------------
    % Outlier Rejection via Correlation
    % ---------------------------------------------------------------------
    % ECG  = ECG - mean(ECG, 2);
    % ECG_avg = mean(ECG, 1);
    %
    % % cdist here represents the Pearson correlation coefficient 'r'
    % cdist = 1 - pdist2(ECG, ECG_avg, 'correlation');
    %
    % % Use strict, absolute thresholds rather than forcing a percentile quota
    % if strcmpi(ECGtype, 'recorded')
    %     threshold_corr = 0.80;
    % else
    %     threshold_corr = 0.70;
    % end
    %
    % % Exclude only the epochs that actually fail the absolute quality check
    % mask_exclude = cdist < threshold_corr;
    % ECGbadEpochs(mask_exclude,:) = [];
    % ECG(mask_exclude,:)          = [];
    % ECGlatency(mask_exclude)     = [];
    %
    % fprintf('Detections removed due to low correlation (threshold %1.2f): %d\n', threshold_corr, sum(mask_exclude));

    % -------------------------------------------------------------------------
    % Find outliers
    % -------------------------------------------------------------------------
    % Generate the time vector centered exactly around 0 (the R-peak)
    T = (0:size(ECG, 2)-1) ./ DATA.srate;
    T = T - mean(T);

    % Find outliers
    [valid_idx, bad_mask, stats] = find_template_outliers(ECG', T);

    ECGEpochs  = ECGEpochs(valid_idx, :);
    ECG        = ECG(valid_idx, :);
    ECGlatency = ECGlatency(valid_idx);

    % Check if any are left
    num_epoch = size(ECG, 1);
    if num_epoch == 0
        warning('Low quality ECG data after artefact rejection... Quitting...');
        ECGmask      = NaN;
        ECGEpochs    = NaN;
        ECGlatency   = NaN;
        heartData    = NaN;
        pulsEstimate = NaN;
        return;
    end

    % ECG mask
    ECGmask = false(size(heartData));
    for i = 1:num_epoch
        ECGmask(ECGEpochs(i, 1):ECGEpochs(i, 2)) = true;
    end

    % ---------------------------------------------------------------------
    % Robust Pulse Estimation
    % ---------------------------------------------------------------------
    IBI = diff(ECGlatency);

    % Filter out physically impossible beat intervals (e.g., under 0.4s or over 2.0s)
    validIBI = IBI(IBI >= 0.4 & IBI <= 2.0);

    % Use median to protect against single missed peaks doubling the interval time
    if ~isempty(validIBI)
        pulsEstimate = median(validIBI) * 60;
    else
        pulsEstimate = NaN;
    end

    fprintf('The final number of heartbeats detected: %d\n', num_epoch);
    fprintf('Average pulse: %1.1f per min.\n', pulsEstimate);

    % ---------------------------------------------------------------------
    % Plotting
    % ---------------------------------------------------------------------
    % Use standard dimensions matching your other preprocessing figures
    fh = figure('Color', 'w', 'Position', [100, 100, 800, 450], 'Visible', optVisible);

    % Explicitly layout using a single compact tile to handle formatting cleanly
    t = tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight');
    ax = nexttile(t);
    hold(ax, 'on');

    % % Generate the time vector centered exactly around 0 (the R-peak)
    % F = (0:size(ECG, 2)-1) ./ DATA.srate;
    % F = F - mean(F);

    % Step 1: Generate a smooth, earthy color palette matrix
    % Using a darker green/teal variant so the faint alpha lines remain visible
    colors = brewermap(num_epoch, 'YlGnBu');

    % Step 2: Plot individual raw epochs with heavy transparency
    % This transforms line clutter into a beautiful statistical density cloud
    for i_ep = 1:num_epoch
        h_line = plot(ax, T, ECG(i_ep, :), 'LineWidth', 0.8, 'Color', colors(i_ep, :));
        h_line.Color(4) = 0.12; % Set alpha channel (opacity) to 12%
    end

    % Step 3: Layer the Grand Average clearly on top
    % A vibrant, solid dark red provides a stark, clear contrast against the background cloud
    ECG_avg = mean(ECG, 1);
    h_avg = plot(ax, T, ECG_avg, 'Color', [0.65, 0.10, 0.10], 'LineWidth', 2.5, 'DisplayName', 'Grand Average');

    % Step 4: Refined Axis Styling & Grid Typography
    grid(ax, 'on');
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Layer', 'top');
    set(ax, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);

    % Scale constraints
    axis(ax, 'tight');
    xlim(ax, [T(1), T(end)]);

    % Clear, non-overlapping labels
    xlabel(ax, 'Time Relative to R-Peak (s)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Amplitude (\muV)', 'FontSize', 12, 'FontWeight', 'bold');

    % Clean Title incorporating metadata dynamically
    title_str = sprintf('Detected ECG QRS (N = %d; %s)', num_epoch, ECGtype);
    title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % Unobtrusive legend pointing only to the average signal trace
    legend(h_avg, 'Grand Average QRS', 'Location', 'NorthEast', 'Box', 'off');

    hold(ax, 'off');

    % Save
    save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_detected_ecg'], [20 11]);
end

end