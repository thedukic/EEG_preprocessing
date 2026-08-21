function [EEG, NumberTrials] = detect_badepochs(EEG, cfg)

% =========================================================================
% Data quaility checks
% =========================================================================

fprintf('\n================================\n');
fprintf('Detecting bad epochs\n');
fprintf('================================\n');

% Note the number of trials
NumberTrials = NaN(4, 1);
NumberTrials(1) = size(EEG.data, 3);

% Find EEG channels
chan_mask = strcmp({EEG.chanlocs.type}, 'EEG');
chan_indx = find(chan_mask);

% =============================================
% A. EEGLAB-based rejection
% =============================================
% Any one of these functions can be commented out to ignore those artifacts when creating the mask
% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs

% fprintf('\n--------------------------------\n');
% fprintf('Max. amplitude (>abs(%d uV))\n', cfg.epoch.amplitude_max);
% fprintf('--------------------------------\n');
% EEG = pop_eegthresh(EEG, 1, ROIidx, -cfg.epoch.amplitude_max, cfg.epoch.amplitude_max, EEG.xmin, EEG.xmax, 1, 0);

fprintf('\n--------------------------------\n');
fprintf('Improbable data\n');
fprintf('--------------------------------\n');
EEG = pop_jointprob(EEG, 1, chan_indx, cfg.epoch.singleChannelImprobableDataThreshold, cfg.epoch.allChannelImprobableDataThreshold, 1, 0);

fprintf('\n--------------------------------\n');
fprintf('Kurtosis\n');
fprintf('--------------------------------\n');
EEG = pop_rejkurt(EEG, 1, chan_indx, cfg.epoch.singleChannelKurtosisThreshold, cfg.epoch.allChannelKurtosisThreshold, 1, 0);

fprintf('\n--------------------------------\n');
fprintf('Max. amplitude (> %d uV)\n', cfg.epoch.amplitude_max);
fprintf('--------------------------------\n');

freq_stop = [2 40]; % Hz
fs = EEG.srate;

% Prototype order 2 -> 4th-order IIR -> 8th-order effective in filtfilt
% Provides strong stopband attenuation without long settling times
[b, a] = butter(2, freq_stop / (fs / 2), 'stop');
assert(isstable(b, a), 'Bandstop filter unstable.');

% Extract channels [channels x timepoints x epochs]
data_eval = double(EEG.data(chan_mask, :, :));
[n_ch, n_pnts, n_trials] = size(data_eval);

% Permute to [timepoints x channels x epochs] for efficient column processing
data_eval = permute(data_eval, [2, 1, 3]);

% 1. Detrend and demean BEFORE filtering to prevent DC step impulse
data_eval = detrend(data_eval);

% 2. Apply zero-phase bandstop filter
data_filt = zeros(size(data_eval));
for i_epoch = 1:n_trials
    data_filt(:, :, i_epoch) = filtfilt(b, a, data_eval(:, :, i_epoch));
end

% 3. Define evaluation window (exclude the outer 50 ms to ignore boundary reflection)
edge_trim_pts = max(1, round(0.050 * fs));
eval_idx = (1 + edge_trim_pts) : (n_pnts - edge_trim_pts);

% 4. Evaluate maximum absolute amplitude across valid time points and channels
% Shape: [channels x epochs]
max_per_chan_trial = squeeze(max(abs(data_filt(eval_idx, :, :)), [], 1));

% 5. Trial rejection mask: flag trial if any channel exceeds threshold
mask_amplitude = any(max_per_chan_trial > cfg.epoch.amplitude_max, 1)';

fprintf('%d/%d trials marked for rejection.\n', sum(mask_amplitude), n_trials);

fprintf('\n--------------------------------\n');
fprintf('EMG envelope\n');
fprintf('--------------------------------\n');

cfg_filt = extract_lpfilt(EEG.ALSUTRECHT.subject.task, cfg);
if isempty(cfg_filt.lp)
    cfg_filt.lp(1) = Inf;
end

filter_emg = 70; % Hz

% Only run EMG envelope detection if data bandwidth permits frequencies >70 Hz
if (cfg_filt.lp(1) - filter_emg) > 20
    data_tmp = double(EEG.data(chan_mask, :, :));
    % Concatenate channels across continuous time if detector expects 2D
    data_tmp_2d = reshape(data_tmp, n_ch, []);

    mask_envelope_pts = detect_emg_envelopes(data_tmp_2d, fs, filter_emg);
    mask_envelope_2d  = reshape(mask_envelope_pts, [n_pnts, n_trials]);
    mask_envelope     = any(mask_envelope_2d, 1)'; % [n_trials x 1] logical

    fprintf('%d/%d trials marked for rejection by EMG envelope.\n', sum(mask_envelope), n_trials);
else
    fprintf('Bypassing EMG envelope check (Data low-pass filtered at %d Hz; requires >%d Hz).\n', ...
        cfg_filt.lp(1), filter_emg + 20);
    mask_envelope = false(n_trials, 1);
end

% EEGTMP = EEG;
% EEGTMP.data(chan_mask, :, mask_envelope) = 0;
% vis_artifacts(EEGTMP, EEG);

% badTrialTmp = find(mask_amplitude);
% figure; tiledlayout(1, 2);
% mytopoplot(mean(data_tmp(:,:,badTrialTmp).^2,[2 3]),[],'Filtered',nexttile); colorbar;
%
% data_tmp = double(EEG.data(1:128,:,:));
% mytopoplot(mean(data_tmp(:,:,badTrialTmp).^2,[2 3]),[],'Raw',nexttile); colorbar;
% disp(sum(mask_amplitude));
%
% [NCHN, NPTS, NTRL]= size(data_tmp);
% for i = 1:NTRL
%     [psdspectra(:, :, i), freq] = pwelch(data_tmp(:,:,i)',NPTS,0,NPTS,EEG.srate);
% end
% figure; plot(freq, mean(psdspectra,3));

fprintf('\n--------------------------------\n');
fprintf('Combining and rejecting\n');
fprintf('--------------------------------\n');

EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
mask_eeglab = EEG.reject.rejglobal;
mask_all = mask_eeglab(:) | mask_amplitude(:) | mask_envelope(:);

EEG = pop_rejepoch(EEG, mask_all, 0);
NumberTrials(2) = EEG.trials;
fprintf('Retained %d/%d trials after Stage A rejection.\n', NumberTrials(2), NumberTrials(1));

% =============================================
% B. EMG-slope-based rejection
% =============================================
fprintf('\n--------------------------------\n');
fprintf('EMG slopes\n');
fprintf('--------------------------------\n');

% Consider targeted MWF?
% EEG = denoise_emg(EEG);

bad_elecpertrial  = cell(1, EEG.trials);
mask_emg          = false(EEG.trials, 1);
report            = struct('listFixed', [], 'listRemove', []);

% % Check the lowpass fitler settings
% % No point of doing this emg check if data are lowpass filtered
% if matches(cfg.epoch.interpolation, {'yes', 'no'}, 'IgnoreCase', true)
%     % Extract lowpass filter settings
%     cfg_filt = extract_lpfilt(EEG.ALSUTRECHT.subject.task, cfg);
%     if ~isempty(cfg_filt.lp)
%         assert(cfg_filt.lp(1) > cfg.ica.emgSlopeFreq2(2), 'Lowpass filtering cutoff too low.');
%     end
%
%     % Find EMG bursts
%     % cfg.bch.emgPeripheral     = false;
%     % cfg.bch.emgSlopeThreshold = 0;
%     % cfg.bch.emgSlopeFreq      = [7 75];
%     % [~, mask_emg] = detect_emg(EEG, cfg.bch);
%
%     chan_mask = strcmp({EEG.chanlocs.type}, 'EEG');
%     chan_locs = EEG.chanlocs(chan_mask);
%
%     [bad_mask, bad_channels, metrics] = detect_emg_unified(EEG.data(chan_mask, : , :), EEG.srate, 'power', cfg);
% end

switch cfg.epoch.interpolation
    case 'yes'
        % % Interpolate
        % fprintf('Interpolating trials with strong EMG leftover.\n');
        %
        % % Inerpolate trials that do not have a lot of electrodes contaminated by EMG
        % if any(bad_mask, "all")
        %     % Organise input for interpolation
        %     bad_elecpertrial = arrayfun(@(x) find(bad_mask(:,x)), 1:EEG.trials, 'UniformOutput', false);
        %
        %     % Max number of contaminated electrodes
        %     fprintf('The maximum number of EMG-contaminated electrodes: %d\n', cfg.epoch.interpolation_max);
        %     [eeg_tmp, report] = interpolate_epochs(EEG.data(chan_mask, :, :), chan_locs, bad_elecpertrial, [], cfg.epoch.interpolation_max);
        %
        %     EEGTMP = EEG;
        %     EEGTMP.data(chan_mask, :, :) = eeg_tmp;
        %     EEGTMP.data(chan_mask, :, report.listRemove) = 0;
        %     vis_artifacts(EEGTMP, EEG);
        %
        %     % Return the data to the struct
        %     EEG.data(chan_mask, :, :) = eeg_tmp;
        %
        %     % These are too noisy to be saved
        %     mask_emg(report.listRemove) = true;
        %
        %     % Check if all trials are still very noisy
        %     if all(mask_emg)
        %         warning('All trials are still full of EMG activity. Consider excluding this participant.'); return;
        %     elseif any(mask_emg)
        %         EEG = pop_rejepoch(EEG, mask_emg, 0);
        %     else
        %         fprintf('Nice, there are no very contaminated trials for rejection.\n');
        %     end
        %
        % else
        %     fprintf('Nice, no EMG found in the data. Skipping trial interpolation.\n');
        % end

    case 'no'
        % % Remove them simply
        % fprintf('Removing trials with strong EMG leftover.\n');
        % mask_emg = any(bad_mask, 1);
        %
        % % Check if all trials are still very noisy
        % if all(mask_emg)
        %     warning('All trials are still full of EMG activity. Consider excluding this participant.'); return;
        % elseif any(mask_emg)
        %     EEG = pop_rejepoch(EEG, mask_emg, 0);
        % else
        %     fprintf('Nice, there are no very contaminated trials for rejection.\n');
        % end

    case 'skip'
        fprintf('Skipping trial-level EMG interpolation/rejection.\n');
end

% Store interpolation metadata
EEG.ALSUTRECHT.epochRejections.InterpTrialInfo = bad_elecpertrial;
EEG.ALSUTRECHT.epochRejections.InterpReport    = report;
EEG.ALSUTRECHT.epochRejections.interpEpochs    = length(report.listFixed);

NumberTrials(3) = EEG.trials;

% =============================================
% C. Detection using variance and the G-ESD method
% =============================================
% fprintf('\n--------------------------------\n');
% fprintf('G-ESD method with variance\n');
% fprintf('--------------------------------\n');
% % EEG = do_gsd(EEG);
% fprintf('Turned off.\n');

% Note the number of trials
NumberTrials(4) = EEG.trials;

% =========================================================================
% Log
EEG.ALSUTRECHT.epochRejections.initialEpochs       = NumberTrials(1);
EEG.ALSUTRECHT.epochRejections.afterEEGLABEpochs   = NumberTrials(2);
EEG.ALSUTRECHT.epochRejections.afterEMGSlopeEpochs = NumberTrials(3);
EEG.ALSUTRECHT.epochRejections.remainingEpochs     = NumberTrials(4);
EEG.ALSUTRECHT.epochRejections.mask_eeglab         = mask_eeglab;
EEG.ALSUTRECHT.epochRejections.mask_amplitude      = mask_amplitude;
EEG.ALSUTRECHT.epochRejections.mask_envelope       = mask_envelope;
EEG.ALSUTRECHT.epochRejections.mask_emg            = mask_emg;
EEG.ALSUTRECHT.epochRejections.proportionOfEpochsRejected = (EEG.ALSUTRECHT.epochRejections.initialEpochs - EEG.ALSUTRECHT.epochRejections.remainingEpochs) / EEG.ALSUTRECHT.epochRejections.initialEpochs;

end

% =========================================================================
% Helper: Extract Task-Specific Low-Pass Filter
% =========================================================================
function cfg_filt = extract_lpfilt(task, cfg)
if strcmpi(task, 'SART') || strcmpi(task, 'MMN')
    cfg_filt = cfg.flt.erp;
elseif strcmpi(task, 'MT')
    cfg_filt = cfg.flt.mt;
elseif any(strcmpi(task, {'RS', 'EO', 'EC'}))
    cfg_filt = cfg.flt.rs;
else
    error('Unknown task specified: %s', task);
end
end