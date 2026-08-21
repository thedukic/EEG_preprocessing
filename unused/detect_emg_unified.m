function [mask, bad_channel, metrics] = detect_emg_unified(data, srate, method, cfg)
% DETECT_EMG_UNIFIED Detects EMG artifacts in continuous or epoched data.
% Inputs:
%   data   - 2D (signals x time) or 3D (signals x time x trials) array
%   srate  - Sampling rate in Hz
%   method - 'slope' or 'asr'
%   cfg    - Configuration structure with threshold parameters

[num_signals, num_pts, num_trials] = size(data);

% -------------------------------------------------------------------------
% 1. Segment Data into 1-second epochs
% -------------------------------------------------------------------------
pts_per_sec = 2*round(srate);
num_segments = floor(num_pts / pts_per_sec);

if num_segments > 1
    data_seg = data(:, 1:(num_segments * pts_per_sec), :);
    data_seg = reshape(data_seg, num_signals, pts_per_sec, num_segments * num_trials);
else
    data_seg = reshape(data, num_signals, num_pts, num_trials);
    num_segments = 1;
end

total_segments = size(data_seg, 3);
mask = false(num_signals, total_segments);
metrics = struct();

% -------------------------------------------------------------------------
% 2. Process based on selected method
% -------------------------------------------------------------------------
switch lower(method)
    case 'asr'
        % % Ensure Nyquist is strictly higher than the maximum defined frequency (55 Hz)
        % nyquist = srate / 2;
        % if nyquist <= 55
        %     error('Sampling rate must be strictly greater than 110 Hz for this ASR filter definition.');
        % end
        %
        % % Both vectors now contain exactly 13 elements (appended 0.1 to m_response)
        % f_response = [0, 2,  3,  7,  8, 12, 14, 20, 25, 45, 50, 55, nyquist];
        % m_response = [0, 0, 0.2, 0.2, 0.5, 0.5, 0.8, 0.8, 1.0, 1.0, 0.1, 0.1, 0.1];
        %
        % try
        %     b_yule = yulewalk(8, f_response / nyquist, m_response);
        %     a_yule = 1;
        % catch ME
        %     % Print the exact internal MATLAB error if it still fails
        %     fprintf('Yulewalk failed: %s\n', ME.message);
        %     [b_yule, a_yule] = butter(4, [15 60]/nyquist, 'bandpass');
        %     warning('Falling back to Butterworth.');
        % end
        %
        % % Plot the filter response
        % figure;
        % freqz(b_yule, a_yule, 1024, srate);
        % title('ASR Heuristic IIR Filter Response');
        %
        % % 2. Apply filter per segment to avoid ringing at spliced boundaries
        % data_filtered = zeros(num_signals, pts_per_sec, total_segments);
        % for i_seg = 1:total_segments
        %     data_filtered(:, :, i_seg) = filtfilt(b_yule, a_yule, data_seg(:, :, i_seg)')';
        % end
        %
        % % 3. Reduce window length to detect transients within the 1-second epoch
        % window_len = round(srate * 0.5);
        %
        % for i_seg = 1:total_segments
        %     for i_sig = 1:num_signals
        %         rms_signal = sqrt(movmean(data_filtered(i_sig, :, i_seg).^2, window_len));
        %
        %         med_rms = median(rms_signal);
        %         mad_rms = median(abs(rms_signal - med_rms)) * 1.4826;
        %         if mad_rms == 0, mad_rms = 1e-6; end
        %
        %         z_profile = (rms_signal - med_rms) / mad_rms;
        %
        %         noisy_indices = z_profile < cfg.emg.asr_z_low | z_profile > cfg.emg.asr_z_high;
        %
        %         if (sum(noisy_indices) / pts_per_sec) > 0.25
        %             mask(i_sig, i_seg) = true;
        %         end
        %     end
        % end

    case 'slope'
        psdspectra = NaN(floor(pts_per_sec/2+1), num_signals, total_segments);
        for i_seg = 1:total_segments
            [psdspectra(:, :, i_seg), freq] = pwelch(data_seg(:, :, i_seg)', pts_per_sec, 0, pts_per_sec, srate);
        end
        psdspectra = permute(psdspectra, [1 3 2]);

        if isscalar(cfg.emg.slope_freq_excl)
            freq_excl = cfg.emg.slope_freq_excl + 2 * [-1 1];
        else
            freq_excl = cfg.emg.slope_freq_excl;
        end

        slope1 = extract_slopes(psdspectra, freq, cfg.emg.slope_freq_1, freq_excl);
        slope2 = extract_slopes(psdspectra, freq, cfg.emg.slope_freq_2, freq_excl);

        % mask = (slope1 > cfg.emg.slope_threshold_1) | (slope2 > cfg.emg.slope_threshold_2);
        % metrics.slopes = max(cat(3, slope1, slope2), [], 3);

        mask = slope1 > cfg.emg.slope_threshold_1;
        metrics.slopes = slope1;

        figure; tiledlayout(2,1);
        nexttile; imagesc(slope1); clim([-5 5]); colorbar;
        nexttile; imagesc(slope2); clim([-8 8]); colorbar;
        figure; imagesc(mask);
    case 'power'

        % =====================================================================
        % Strong EMG
        % =====================================================================
        % 2. EEG: Temporarily highpass filter
        [bh, ah] = butter(4, 70/(srate/2), 'high');
        dataeeg = data_seg(:, :);
        dataeeg = filtfilt(bh, ah, dataeeg')';
        % for i_seg = 1:total_segments
        %     dataeeg(:, :, i_seg) = filtfilt(b_yule, a_yule, data_seg(:, :, i_seg)')';
        % end

        % Traces of EMG power
        dataeeg = 15 * abs(dataeeg);

        % ======================================
        % Try to find clusters of higher numbers
        % ======================================
        P1 = 99;
        treshold = prctile(dataeeg(:), P1);
        maskEMG = dataeeg > treshold;

        data1 = dataeeg;
        data1(~maskEMG) = 0;
        data1(maskEMG)  = 1;

        smoothLengthSample = @(smoothLengthMS) round(srate * (smoothLengthMS/1000));
        P = smoothLengthSample(500);
        data1 = movmean(data1', P)';

        % ======================================
        % At least X% of electrodes must be affected
        % ======================================
        P2 = 95;
        dataTmp = data1;
        dataTmp(dataTmp == 0) = [];
        treshold = prctile(dataTmp, P2);

        mask = data1 > treshold;
        mask = 100 * mean(mask,1);

        % Mask1
        % -> Many channels affected together with the EOG
        % -> This happens due in large EMG/movement/blink artifacts
        P3 = 15; % X% of elec
        extremeMaskTmp2 = mask >= 2*P3;

        % Mask
        mask(mask < P3) = 0;
        mask(extremeMaskTmp2) = 2*P3;

        P = smoothLengthSample(1000);
        mask = movmean(mask, P) > 0;
        mask = reshape(mask, [], total_segments);

end

% -------------------------------------------------------------------------
% 3. Contamination Logic & Reconstruct Original Dimensions
% -------------------------------------------------------------------------
if num_segments > 1
    segment_mask_3d = reshape(mask, num_signals, num_segments, num_trials);

    % Use reshape to prevent dimension collapse if N=1 or num_signals=1
    bad_ratio = squeeze(mean(segment_mask_3d, 2));
    mask = bad_ratio > (1/3);
end

% Average time is greater than set
bad_channel = mean(mask, 2) > cfg.emg.slope_time;

% indx = find(any(mask, 1));
% figure; tiledlayout(1, 2);
% for i = 1:length(indx)
%     mask_tmp = mask(:, indx(i));
%     data_tmp = squeeze(psdspectra(:, indx(i), mask_tmp));
%     nexttile; plot(log10(freq), log10(data_tmp)); axis tight;
%     data_tmp = squeeze(psdspectra(:, indx(i) + 1, mask_tmp));
%     nexttile; plot(log10(freq), log10(data_tmp)); axis tight;
%     pause; clf;
% end

end

% =========================================================================
% HELPER: Vectorised Slope Extraction
% =========================================================================
function slopes = extract_slopes(psdspectra, freq_axis, freq_range, freq_excl)
[num_freqs, num_trials, num_signals] = size(psdspectra);

mask_include = freq_axis >= freq_range(1) & freq_axis <= freq_range(2);
mask_exclude = freq_axis >= freq_excl(1) & freq_axis <= freq_excl(2);
mask_include = mask_include & ~mask_exclude;

logpow = log10(psdspectra(mask_include, :, :));
logfoi = log10(freq_axis(mask_include));

X = [logfoi(:), ones(length(logfoi), 1)];
logpow_2d = reshape(logpow, size(logpow, 1), []);
P = X \ logpow_2d;

slopes = reshape(P(1, :), num_trials, num_signals)';
end