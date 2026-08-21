function EEG = estimate_channelcov(EEG, cfg)
% ESTIMATE_CHANNELCOV
% Vectorised computation of broadband and high-frequency (55-95 Hz) channel
% covariance matrices across scalp leads and bipolar EOGs.
%
% Inputs:
%   EEG : Cleaned EEGLAB dataset
%   cfg : (Optional) Configuration structure:
%         cfg.hf_band : High-frequency band in Hz [Default: [55, 95]]
%
% Outputs:
%   EEG : Updated dataset with matrices and scalar QC metrics logged under
%         EEG.ALSUTRECHT.channelcov

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'hf_band'), cfg.hf_band = [55, 95]; end

fs = EEG.srate;

% -------------------------------------------------------------------------
% 1. Select Channels (128 Scalp EEG + 2 EOG = 130 Channels)
% -------------------------------------------------------------------------
eegchan = strcmp({EEG.chanlocs.type}, 'EEG');
eogchan = ismember({EEG.chanlocs.labels}, {'VEOG', 'HEOG'});
selchan = eegchan | eogchan;

num_channels = sum(selchan);
assert(num_channels == 130, 'Expected 130 channels, found %d', num_channels);

% Boolean mask for scalp channels within the selected subset
scalp_mask = eegchan(selchan);

% Extract selected channels in double precision: [channels x pnts x trials]
D_3d = double(EEG.data(selchan, :, :));
[n_ch, n_pnts, n_trials] = size(D_3d);

% -------------------------------------------------------------------------
% 2. Broadband Covariance Matrix
% -------------------------------------------------------------------------
% Flatten across time and trials -> [Total Samples x Channels]
D_raw_2d  = reshape(D_3d, n_ch, [])';
C_cov_raw = cov(D_raw_2d, 'omitrows');

% -------------------------------------------------------------------------
% 3. High-Frequency (55-95 Hz) Covariance (Filtered Per Epoch)
% -------------------------------------------------------------------------
% 4th-order zero-phase Butterworth bandpass filter
[b_hf, a_hf] = butter(4, cfg.hf_band / (fs / 2), 'bandpass');

% Permute to [pnts x channels x trials] so filtfilt operates down time (dim 1)
% Filtering per epoch avoids boundary step ringing across concatenated trials
D_perm    = permute(D_3d, [2, 1, 3]);
D_hf_perm = zeros(size(D_perm));

for tr = 1:n_trials
    tr_data = D_perm(:, :, tr);
    if ~any(isnan(tr_data(:)))
        % Filters all 130 channels simultaneously for this epoch
        D_hf_perm(:, :, tr) = filtfilt(b_hf, a_hf, tr_data);
    else
        D_hf_perm(:, :, tr) = NaN;
    end
end

% Reshape filtered data -> [Total Samples x Channels]
D_hf_2d  = reshape(permute(D_hf_perm, [2, 1, 3]), n_ch, [])';
C_cov_hf = cov(D_hf_2d, 'omitrows');

% -------------------------------------------------------------------------
% 4. Scalar Quality Control Metrics for Group Outlier Detection
% -------------------------------------------------------------------------
var_raw = diag(C_cov_raw);
var_hf  = diag(C_cov_hf);

% Isolate scalp-only variances
var_raw_scalp = var_raw(scalp_mask);
var_hf_scalp  = var_hf(scalp_mask);

% Compute temporary correlation on the fly to get inter-channel coupling metric
C_corr_hf     = corrcov(C_cov_hf);
off_diag_mask = ~eye(num_channels);

qc.total_broadband_power = sum(var_raw_scalp);
qc.total_hf_power        = sum(var_hf_scalp);
qc.mean_hf_variance      = mean(var_hf_scalp);
qc.max_hf_channel_ratio  = max(var_hf_scalp) / median(var_hf_scalp);
qc.mean_hf_correlation   = mean(abs(C_corr_hf(off_diag_mask)));

% -------------------------------------------------------------------------
% 5. Log Output to EEG Structure
% -------------------------------------------------------------------------
EEG.ALSUTRECHT.channelcov.cov_clean  = C_cov_raw;
EEG.ALSUTRECHT.channelcov.cov_hf     = C_cov_hf;
EEG.ALSUTRECHT.channelcov.qc_metrics = qc;

end