function EEG = detect_ic_asr_transients(EEG, data, cfg)
% DETECT_IC_ASR_TRANSIENTS Identifies ICs dominated by high-amplitude 
% transient artifacts using an adapted ASR window-rejection heuristic.

fprintf('\n--------------------------------\n');
fprintf('ASR Transient IC Detection\n');
fprintf('--------------------------------\n');

% Set up default parameters if not defined in cfg
if ~isfield(cfg.thresh, 'asr_z_low'),  cfg.thresh.asr_z_low = -3.5; end
if ~isfield(cfg.thresh, 'asr_z_high'), cfg.thresh.asr_z_high = 7.0; end

% -------------------------------------------------------------------------
% 1. Apply the Custom Heuristic IIR Spectral Weighting
% -------------------------------------------------------------------------
% The Yule-Walker filter emphasizes frequencies where muscle/transient noise 
% dominates while suppressing clean alpha/EEG rhythms.
% Designing an 8th-order Yule-Walker filter based on the ASR architecture:
f_response = [0, 2,  3,  7,  8, 12, 14, 20, 25, 45, 50, 55, size(data, 2)]; 
m_response = [0, 0, 0.2, 0.2, 0.5, 0.5, 0.8, 0.8, 1.0, 1.0, 0.1, 0.1]; % Suppression around line noise

% Safe fallback to a standard high-order bandpass if yulewalk parameters vary
try
    b_yule = yulewalk(8, f_response / (EEG.srate/2), m_response);
    a_yule = 1; % FIR approximation or design dependent
    data_filtered = filtfilt(b_yule, a_yule, data')';
catch
    % Clean fallback to the main transient passband (15 - 60 Hz)
    [b_yule, a_yule] = butter(4, [15 60]/(EEG.srate/2), 'bandpass');
    data_filtered = filtfilt(b_yule, a_yule, data')';
end

% -------------------------------------------------------------------------
% 2. Calculate Moving RMS Amplitude & 3. Perform Z-Score Transform
% -------------------------------------------------------------------------
% Using a standard 500 ms moving window to match typical ASR calibration scales
window_len = round(EEG.srate * 1); 
num_samples = size(data_filtered, 2);

% Preallocate flagging arrays
num_ics = size(data, 1);
bad_ic_mask = false(num_ics, 1);
outlier_ratio = zeros(num_ics, 1);

for i_ic = 1:num_ics
    % Root-Mean-Square (RMS) amplitude extraction via a moving vector
    ic_data = data_filtered(i_ic, :);
    rms_signal = sqrt(movmean(ic_data.^2, window_len));
    
    % Robust Z-score transform based on median/MAD to prevent the massive 
    % transient spikes from skewing the mean calculation itself
    med_rms = median(rms_signal);
    mad_rms = median(abs(rms_signal - med_rms)) * 1.4826;
    if mad_rms == 0, mad_rms = 1e-6; end
    
    z_profile = (rms_signal - med_rms) / mad_rms;
    
    % ---------------------------------------------------------------------
    % 4. Find Datapoints Outside -3.5 and 7.0
    % ---------------------------------------------------------------------
    noisy_indices = z_profile < cfg.thresh.asr_z_low | z_profile > cfg.thresh.asr_z_high;
    
    % If more than 3% of the total recording timeline is dominated by these 
    % extreme asymmetric spikes, the component is deemed non-stationary.
    outlier_ratio(i_ic) = sum(noisy_indices) / num_samples;
    if outlier_ratio(i_ic) > 0.03
        bad_ic_mask(i_ic) = true;
    end
end


end