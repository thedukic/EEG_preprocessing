function [emg_ratio, fail_spectral, chan_ratios] = estimate_noise_ratio(EEG, f_low, f_high, thresh)
% ESTIMATE_NOISE_RATIO Evaluates high-frequency muscle noise relative 
% to low-frequency neural rhythmicity on epoched EEG data.
%
% Syntax:
%   [emg_ratio, fail_spectral, chan_ratios] = estimate_noise_ratio(EEG, f_low, f_high, thresh)
%
% Inputs:
%   EEG           - Epoched EEGLAB dataset structure [channels x time x trials]
%   f_low         - [1x2] Low-frequency neural band in Hz (default: [4 12])
%   f_high        - [1x2] High-frequency muscle band in Hz (default: [30 45])
%   thresh        - Ratio threshold for flagging severe contamination (default: 0.40)
%
% Outputs:
%   emg_ratio     - Whole-scalp mean ratio (High / Low power)
%   fail_spectral - Logical flag (true if emg_ratio > thresh)
%   chan_ratios   - [1 x n_chans] Ratio per channel to evaluate spatial focus

% -------------------------------------------------------------------------
% 1. Input Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(f_low),  f_low  = [4 12];  end
if nargin < 3 || isempty(f_high), f_high = [30 45]; end
if nargin < 4 || isempty(thresh), thresh = 0.40;    end

assert(isstruct(EEG) && isfield(EEG, 'data'), 'Input must be a valid EEGLAB structure.');
assert(ndims(EEG.data) == 3, 'Input EEG.data must be 3D epoched data [chans x time x trials].');

% Identify EEG channels safely
if isfield(EEG.chanlocs, 'type') && any(strcmpi({EEG.chanlocs.type}, 'EEG'))
    mask_eeg = strcmpi({EEG.chanlocs.type}, 'EEG');
else
    mask_eeg = true(1, EEG.nbchan);
end

% -------------------------------------------------------------------------
% 2. Power Spectrum Estimation (Trial-by-Trial)
% -------------------------------------------------------------------------
eeg_data = double(EEG.data(mask_eeg, :, :));
[n_chans, n_pnts, n_trials] = size(eeg_data);

% Demean each trial along the time axis
eeg_data = eeg_data - mean(eeg_data, 2);

fs = EEG.srate;
win_len = min(n_pnts, 2 * fs); % 2-second window or full epoch length
noverlap = 0;
nfft = win_len;

% Initialise with trial 1 to obtain frequency bins
[pxx_init, f] = pwelch(eeg_data(:, :, 1)', win_len, noverlap, nfft, fs);
psd_trials = zeros(length(f), n_chans, n_trials);
psd_trials(:, :, 1) = pxx_init;

for i_trl = 2:n_trials
    psd_trials(:, :, i_trl) = pwelch(eeg_data(:, :, i_trl)', win_len, noverlap, nfft, fs);
end

% Mean PSD across all clean trials [frequencies x channels]
pxx = mean(psd_trials, 3, 'omitnan');

% -------------------------------------------------------------------------
% 3. Band Power Integration and Ratio Calculation
% -------------------------------------------------------------------------
idx_low  = f >= f_low(1)  & f <= f_low(2);
idx_high = f >= f_high(1) & f <= f_high(2);

% Integrate power across bands (\muV^2)
p_low_chan  = trapz(f(idx_low),  pxx(idx_low, :), 1);
p_high_chan = trapz(f(idx_high), pxx(idx_high, :), 1);

% Compute ratio per channel and whole-scalp mean
chan_ratios   = p_high_chan ./ (p_low_chan + eps);
emg_ratio     = mean(chan_ratios, 'omitnan');
fail_spectral = emg_ratio > thresh;

% -------------------------------------------------------------------------
% 4. Console Audit Output
% -------------------------------------------------------------------------
status_str = 'PASS';
if fail_spectral
    status_str = 'FLAG (Excessive HF Noise)';
end

fprintf('Spectral Noise Audit [%d-%d Hz / %d-%d Hz]:\n', ...
    f_high(1), f_high(2), f_low(1), f_low(2));
fprintf('  Mean Ratio = %.3f (Thresh: %.2f) | Max Chan = %.3f -> %s\n', ...
    emg_ratio, thresh, max(chan_ratios, [], 'omitnan'), status_str);

end