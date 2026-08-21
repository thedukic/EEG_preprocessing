function [slopes, mask, other] = detect_emg(EEG, cfg)
% =========================================================================
% DETECT_EMG: Sub-segment spectral slope evaluation for EMG detection
% =========================================================================
% Evaluates high-frequency spectral flattening across sub-segments within
% each trial to detect transient muscle bursts in long trials (e.g. 15s MT).
% =========================================================================

fprintf('Detecting EMG in EEG data...\n');

% 1. Extract EEG channels
chaneeg = strcmp({EEG.chanlocs.type}, 'EEG');
num_channels_eeg = sum(chaneeg);
assert(num_channels_eeg > 0, 'No EEG channels found with type == ''EEG''.');

fs = EEG.srate;

% 2. Sub-segment duration configuration (default: 1.0 second windows)
if isfield(cfg.emg, 'seg_len_sec') && ~isempty(cfg.emg.seg_len_sec)
    seg_dur = cfg.emg.seg_len_sec;
else
    seg_dur = 1.0; % 1-second sub-segments
end
seg_pts = round(seg_dur * fs);

% -------------------------------------------------------------------------
% 3. Windowing and Reshaping for Long (15s) vs Short Epochs
% -------------------------------------------------------------------------
if ndims(EEG.data) == 3
    L = EEG.pnts;
    N = EEG.trials;

    num_segments = floor(L / seg_pts);

    if num_segments >= 1
        % Truncate to exact integer multiple of sub-segment length
        valid_pnts = num_segments * seg_pts;
        dataeeg = double(EEG.data(chaneeg, 1:valid_pnts, :));

        % Reshape: [Channels x SegmentPoints x Segments x Trials]
        dataeeg = reshape(dataeeg, num_channels_eeg, seg_pts, num_segments, N);

        % Flatten to [Channels x SegmentPoints x (Segments * Trials)]
        NTRL_processed = num_segments * N;
        dataeeg = reshape(dataeeg, num_channels_eeg, seg_pts, NTRL_processed);
        NPTS = seg_pts;
    else
        % Trial is shorter than seg_pts: evaluate the whole epoch as 1 segment
        dataeeg = double(EEG.data(chaneeg, :, :));
        NPTS = L;
        NTRL_processed = N;
        num_segments = 1;
    end
else
    % Continuous 2D data
    assert(ismatrix(EEG.data), 'EEG.data must be 2D or 3D.');
    data_raw = double(EEG.data(chaneeg, :));

    NPTS = seg_pts;
    N = floor(size(data_raw, 2) / NPTS);
    dataeeg = reshape(data_raw(:, 1:(N * NPTS)), num_channels_eeg, NPTS, N);

    num_segments = 1;
    NTRL_processed = N;
    L = NPTS;
end

other.modulus = size(EEG.data(:, :), 2) - (N * L);
other.L = L;
other.N = N;
other.num_segments = num_segments;

fprintf('Processing %d sub-segments (%d trials x %d sub-segments of %.1fs)\n', ...
    NTRL_processed, N, num_segments, NPTS / fs);

% -------------------------------------------------------------------------
% 4. Fast Vectorised Spectral Estimation Across All Sub-Segments
% -------------------------------------------------------------------------
% Detrend and apply Hanning window across time dimension
w = hanning(NPTS);
w = w / sqrt(mean(w.^2)); % Normalise window energy

% Permute to [TimePoints x (Channels * SubSegments)]
data_2d = reshape(permute(dataeeg, [2, 1, 3]), NPTS, []);
data_detrended = detrend(data_2d);
data_windowed  = data_detrended .* w;

% Single-sided power spectrum
n_fft = 2^nextpow2(NPTS);
fft_out = fft(data_windowed, n_fft, 1);
fft_half = fft_out(1:(n_fft/2 + 1), :);

psdspectra = (abs(fft_half).^2) / (NPTS * fs);
psdspectra(2:end-1, :) = psdspectra(2:end-1, :) * 2;

freq = (0:(n_fft/2)) * (fs / n_fft);

% -------------------------------------------------------------------------
% 5. Linear Regression in Log-Log Space
% -------------------------------------------------------------------------
% Default frequency fitting band (20 to 45 Hz keeps inside 60 Hz LP passband)
if ~isfield(cfg.emg, 'slope_freq_1') || isempty(cfg.emg.slope_freq_1)
    cfg.emg.slope_freq_1 = [20 45];
end

if ~isfield(cfg.emg, 'slope_threshold_1') || isempty(cfg.emg.slope_threshold_1)
    cfg.emg.slope_threshold_1 = -0.50;
end

if ~isfield(cfg.emg, 'segment_ratio_threshold') || isempty(cfg.emg.segment_ratio_threshold)
    cfg.emg.segment_ratio_threshold = 1/3; % Contaminated if > 33% of sub-segments are bad
end

frqmsk = (freq >= cfg.emg.slope_freq_1(1)) & (freq <= cfg.emg.slope_freq_1(2));
assert(sum(frqmsk) >= 3, 'Insufficient frequency points for regression in the selected band.');

logfoi = log10(freq(frqmsk))';
logpow = log10(psdspectra(frqmsk, :) + eps);

% Linear fit: log10(Power) = slope * log10(Freq) + intercept
X = [logfoi, ones(length(logfoi), 1)];
P = X \ logpow;

slopes_vector = P(1, :);
slopesChannelsxProcessed = reshape(slopes_vector, num_channels_eeg, NTRL_processed);

% -------------------------------------------------------------------------
% 6. Aggregate Sub-Segments Back to Original Trial Dimensions
% -------------------------------------------------------------------------
if num_segments > 1
    % Reshape back to [Channels x SubSegments x Original Trials]
    slopes = reshape(slopesChannelsxProcessed, num_channels_eeg, num_segments, N);

    % Identify which sub-segments exceed the slope threshold
    segment_contam = slopes > cfg.emg.slope_threshold_1;

    % Proportion of contaminated sub-segments per channel per trial [Channels x Trials]
    prop_bad_segments = squeeze(mean(segment_contam, 2));

    % Flag channel-trial pair as bad if contaminated proportion exceeds threshold (e.g. > 1/3)
    mask = prop_bad_segments > cfg.emg.segment_ratio_threshold;

    other.prop_bad_segments = prop_bad_segments;
else
    slopes = slopesChannelsxProcessed;
    mask = slopes > cfg.emg.slope_threshold_1;
    other.prop_bad_segments = double(mask);
end

end





% function [slopes, other] = detect_emg(EEG, cfg)
% % The original reference uses 7-75 Hz and the RELAX toolbox as well, slope treshold >-0.31 or -0.51
% % But MNE toolbox suggests 7-45 Hz and says that this is a better value in practice, but slope treshold is >0
% % https://mne.tools/dev/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.find_bads_muscle
%
% % Looking into EEG data only
% fprintf('Detecting EMG in EEG data...\n');
% chaneeg = strcmp({EEG.chanlocs.type}, 'EEG');
%
% if ndims(EEG.data) == 3
%     dataeeg = EEG.data(chaneeg, :, :);
%     L = EEG.pnts;
%     N = EEG.trials;
% else
%     dataeeg = EEG.data(chaneeg, :);
%
%     % Epoch into 1s
%     L = EEG.srate;
%     N = floor(size(dataeeg, 2) / L);
%     dataeeg = reshape(dataeeg(:, 1:N*L), num_channels_eeg, L, N);
% end
%
% % Check
% modulus = size(EEG.data(:, :), 2) - N*L;
% assert(modulus >= 0);
%
% % Compute (log) power spectra
% [NCHN, NPTS, NTRL] = size(dataeeg);
% psdspectra = NaN(floor(NPTS/2+1), NCHN, NTRL);
%
% % Fixed the loop variable name to match the iteration limit
% for i_trial = 1:NTRL
%     [psdspectra(:, :, i_trial), freq] = pwelch(dataeeg(:, :, i_trial)', NPTS, 0, NPTS, EEG.srate);
% end
%
% fprintf('Slope frequency range: %d-%d Hz\n', cfg.emgSlopeFreq);
%
% % Transpose to [Frequencies x Trials x Channels]
% psdspectra = permute(psdspectra, [1 3 2]);
%
% frqmsk = freq >= cfg.emgSlopeFreq(1) & freq <= cfg.emgSlopeFreq(2);
% logpow = log10(psdspectra(frqmsk, :, :));
% logfoi = log10(freq(frqmsk));
%
% % -------------------------------------------------------------------------
% % Vectorised linear regression (replaces the nested polyfit loops)
% % -------------------------------------------------------------------------
% % Create the design matrix for a 1st-degree polynomial fit
% X = [logfoi(:), ones(length(logfoi), 1)];
%
% % Flatten logpow to a 2D matrix: [Frequencies x (Trials * Channels)]
% logpow_2d = reshape(logpow, size(logpow, 1), []);
%
% % Solve for all slopes simultaneously using matrix left division
% P = X \ logpow_2d;
%
% % Extract the slopes (first row) and reshape back to [Channels x Epochs]
% slopes_vector = P(1, :);
% slopes = reshape(slopes_vector, NTRL, NCHN)';
% % -------------------------------------------------------------------------
%
% other.modulus = modulus;
% other.L = L;
% other.N = N;
%
% % Exclude electrodes that are unlikely contaminated by EMG
% if cfg.emgPeripheral
%     is_perielec = select_peripheralelecs(EEG);
%
%     % Prevent matrix dimension crash by matching the mask to the EEG subset
%     is_perielec_eeg = is_perielec(chaneeg);
%
%     slopes(~is_perielec_eeg, :) = NaN;
%     fprintf('Using only the slopes of the peripheral electrodes (N = %d).\n', sum(is_perielec_eeg));
% end
%
% end