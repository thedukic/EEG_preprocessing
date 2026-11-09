function EEG = filter_signal(EEG, lp, hp, chans, data_type)
% FILTER_SIGNAL: Filters EEG data (EEGLAB or FieldTrip format) with mirror padding.
% Arguments:
%    EEG       - EEG data structure (EEGLAB/FieldTrip format).
%    lp        - Low-pass filter parameters: [Cutoff (Hz), Order].
%    hp        - High-pass filter parameters: [Cutoff (Hz), Order].
%    chans     - Channels to filter (indices).
%    pad       - Duration of mirror padding in seconds (default: 2.0 s).
%    data_type - Data type: 'eeglab' or 'fieldtrip'.

% =========================================================================
% !!! IMPORTANT !!!
% Because filtfilt performs a forward-reverse filtering, it doubles the filter order!
% Filters should be applied to continuous (or boundary-delimited) data chunks.
%
% Order 1 ~ 6 dB/oct  -> becomes Order 2 (12 dB/oct)
% Order 2 ~ 12 dB/oct -> becomes Order 4 (24 dB/oct)
% Order 4 ~ 24 dB/oct -> becomes Order 8 (48 dB/oct)
% =========================================================================

% Check Inputs
pad     = 1.0; % Default X seconds of mirror reflection
dur_min = 0.9; % Minimum duration of data that will be filtered - mayeb this should be zero?

if ~isempty(lp)
    validateattributes(lp, {'numeric'}, {'numel', 2, 'positive'}, mfilename, 'lp');
else
    lp = [];
end
if ~isempty(hp)
    validateattributes(hp, {'numeric'}, {'numel', 2, 'positive'}, mfilename, 'hp');
else
    hp = [];
end
assert(ismember(lower(data_type), {'fieldtrip', 'eeglab'}), 'Unsupported data type.');

% Define variables and Nyquist frequency
if strcmpi(data_type, 'fieldtrip')
    num_trial   = length(EEG.trial);
    srate       = EEG.fsample;
    freq_nyq    = srate / 2;
    num_channel = size(EEG.trial{1}, 1);
    warning('Filtering ALL given signals across %d trials...\n', num_trial);
elseif strcmpi(data_type, 'eeglab')
    num_block   = length(EEG);
    srate       = EEG(1).srate;
    freq_nyq    = srate / 2;
    num_channel = length(chans);
    fprintf('Detected: %d block(s)...\n', num_block);
end

% Calculate filter coefficients
filter_type = [~isempty(lp), ~isempty(hp)];
if ~any(filter_type)
    fprintf('Skipping filtering as parameters are not defined.\n');
    return;
end

% Filter design
bLow = []; aLow = [];
bHigh = []; aHigh = [];

if filter_type(1)
    [bLow, aLow] = butter(lp(2)/2, lp(1)/freq_nyq, 'low');
    assert(isstable(bLow, aLow), 'Low-pass filter unstable.');
    fprintf('Using Butterworth lowpass [%1.2f Hz, order %d] filter.\n', lp(1), lp(2));
    if lp(1) < 20, warning('N.B. Very low (<20 Hz) lowpass filter settings.'); end
end

if filter_type(2)
    [bHigh, aHigh] = butter(hp(2)/2, hp(1)/freq_nyq, 'high');
    assert(isstable(bHigh, aHigh), 'High-pass filter unstable.');
    fprintf('Using Butterworth highpass [%1.2f Hz, order %d] filter.\n', hp(1), hp(2));
    if hp(1) > 1.5, warning('N.B. Very high (>1.5 Hz) highpass filter settings.'); end
end

% Padding length in samples
nPadDefault = round(pad * srate);
nDurDefault = round(dur_min * srate);

% Execution: FieldTrip
if strcmpi(data_type, 'fieldtrip')
    for i = 1:num_trial
        assert(size(EEG.trial{i}, 1) == num_channel);

        data_tmp = double(EEG.trial{i});
        EEG.trial{i} = apply_filter_with_mirror_padding(data_tmp, bLow, aLow, bHigh, aHigh, filter_type, freq_nyq, nPadDefault);

        assert(size(EEG.trial{i}, 1) == num_channel);
    end

    % Execution: EEGLAB
else
    for i = 1:num_block
        assert(size(EEG(i).data(chans, :), 1) == num_channel);

        if ~isempty(EEG(i).event)
            eventLabels = {EEG(i).event(:).type};
            eventTimes  = [EEG(i).event(:).latency];
            mask        = contains(eventLabels, 'boundary');
            boundaries  = unique(floor([0, eventTimes(mask), size(EEG(i).data, 2)]));
        else
            boundaries  = [0, size(EEG(i).data, 2)];
        end

        num_segment = length(boundaries) - 1;
        if num_segment > 1
            fprintf('Block %d: Boundaries detected! Each chunk (N = %d) will be filtered separately.\n', i, num_segment);
        else
            fprintf('Block %d: Boundaries not detected.\n', i);
        end

        for j = 1:num_segment
            data_indx = (boundaries(j) + 1):boundaries(j + 1);

            if length(data_indx) > nDurDefault
                data_tmp = double(EEG(i).data(chans, data_indx));

                % Filter with boundary-safe mirror reflection
                data_tmp = apply_filter_with_mirror_padding(data_tmp, bLow, aLow, bHigh, aHigh, filter_type, freq_nyq, nPadDefault);

                % Place back
                EEG(i).data(chans, data_indx) = data_tmp;
            else
                warning('Chunk %d: Too small (<= %d) for filtering (L = %d samples)!', j, nDurDefault, length(data_indx));
            end
        end

        assert(size(EEG(i).data(chans, :), 1) == num_channel);
    end
end

fprintf('Done!\n');

end


% =========================================================================
% HELPER: MIRROR PADDING, DC DETREND, AND ZERO-PHASE FILTERING
% =========================================================================
function data_out = apply_filter_with_mirror_padding(data_in, bLow, aLow, bHigh, aHigh, typeFilter, freq_nyq, nPadReq)
% dataIn: [n_channels x n_samples]

num_sample = size(data_in, 2);

% Ensure padding length does not exceed available chunk length
num_pad = min(nPadReq, num_sample - 1);

% 1. Remove DC offset first so edge mirroring does not amplify baseline steps
data_demean = remove_dcsignal(data_in, freq_nyq);

% 2. Apply mirror (reflection) padding symmetrically
if num_pad > 0
    pad_left    = fliplr(data_demean(:, 1:num_pad));
    pad_right   = fliplr(data_demean(:, (end - num_pad + 1):end));
    data_padded = [pad_left, data_demean, pad_right];
else
    data_padded = data_demean;
end

% 3. Zero-phase filtering on the padded segment
% Transpose so filtfilt operates along the time dimension (columns)
if all(typeFilter)
    data_padded = filtfilt(bLow, aLow, data_padded')';
    data_padded = filtfilt(bHigh, aHigh, data_padded')';
elseif typeFilter(1)
    data_padded = filtfilt(bLow, aLow, data_padded')';
elseif typeFilter(2)
    data_padded = filtfilt(bHigh, aHigh, data_padded')';
end

% 4. Trim off the mirror padding
if num_pad > 0
    data_out = data_padded(:, (num_pad + 1):(end - num_pad));
else
    data_out = data_padded;
end

end