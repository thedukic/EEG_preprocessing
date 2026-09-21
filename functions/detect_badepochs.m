function [EEG, num_trials] = detect_badepochs(EEG, tag, cfg)

% =========================================================================
% Data quaility checks
% =========================================================================

fprintf('\n================================\n');
fprintf('Detecting bad epochs\n');
fprintf('================================\n');

% Note the number of trials
num_trials = NaN(4, 1);
num_trials(1) = size(EEG.data, 3);

% =============================================
% A. Rejection
% =============================================
% Any one of these functions can be commented out to ignore those artifacts when creating the mask
% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs
% EEGLAB detection - I would NOT recommand as it flags true alpha/beta oscillations as bad data!

% % Find EEG channels
% chan_mask = strcmp({EEG.chanlocs.type}, 'EEG');
% chan_indx = find(chan_mask);

% fprintf('\n--------------------------------\n');
% fprintf('Max. amplitude (>abs(%d uV))\n', cfg.epoch.amplitude_max);
% fprintf('--------------------------------\n');
% EEG = pop_eegthresh(EEG, 1, ROIidx, -cfg.epoch.amplitude_max, cfg.epoch.amplitude_max, EEG.xmin, EEG.xmax, 1, 0);

% fprintf('\n--------------------------------\n');
% fprintf('Improbable data\n');
% fprintf('--------------------------------\n');
% EEG = pop_jointprob(EEG, 1, chan_indx, cfg.epoch.singleChannelImprobableDataThreshold, cfg.epoch.allChannelImprobableDataThreshold, 1, 0);
%
% fprintf('\n--------------------------------\n');
% fprintf('Kurtosis\n');
% fprintf('--------------------------------\n');
% EEG = pop_rejkurt(EEG, 1, chan_indx, cfg.epoch.singleChannelKurtosisThreshold, cfg.epoch.allChannelKurtosisThreshold, 1, 0);

fprintf('\n--------------------------------\n');
fprintf('Max. amplitude (> %d uV)\n', cfg.epoch.amplitude_max);
fprintf('--------------------------------\n');
mask_amplitude = detect_high_amplitude(EEG, cfg);

fprintf('\n--------------------------------\n');
fprintf('EMG contamination\n');
fprintf('--------------------------------\n');
mask_muscle_env = detect_high_muscle(EEG, cfg);

fprintf('\n--------------------------------\n');
fprintf('Combining and rejecting\n');
fprintf('--------------------------------\n');
% % Combine EEGLAB
% EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
% mask_eeglab = EEG.reject.rejglobal;
mask_eeglab = false(size(mask_amplitude));

% Combine all
mask_all = mask_eeglab(:) | mask_amplitude(:) | mask_muscle_env(:);

% Checks
if any(mask_all)
    [reject_stats, fh] = check_rejected_trials(EEG, find(mask_all));
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, sprintf('%s_power_removed_epochs_%s', EEG.ALSUTRECHT.subject.id, num2str(tag)), [32 15]);
else
    reject_stats = [];
end

% Remove
EEG = pop_rejepoch(EEG, mask_all, 0);

num_trials(2) = EEG.trials;
fprintf('Retained %d/%d trials after Stage A rejection.\n', num_trials(2), num_trials(1));

% =============================================
% B. Interpolation (or rejection)
% =============================================
fprintf('\n--------------------------------\n');
fprintf('EMG slopes\n');
fprintf('--------------------------------\n');

% Consider targeted MWF?
% EEG = denoise_emg(EEG);

bad_elec_pertrial  = cell(1, EEG.trials);
mask_muscle_slope = false(EEG.trials, 1);
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
        %     bad_elec_pertrial = arrayfun(@(x) find(bad_mask(:,x)), 1:EEG.trials, 'UniformOutput', false);
        %
        %     % Max number of contaminated electrodes
        %     fprintf('The maximum number of EMG-contaminated electrodes: %d\n', cfg.epoch.interpolation_max);
        %     [eeg_tmp, report] = interpolate_epochs(EEG.data(chan_mask, :, :), chan_locs, bad_elec_pertrial, [], cfg.epoch.interpolation_max);
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
EEG.ALSUTRECHT.epochRejections.InterpTrialInfo = bad_elec_pertrial;
EEG.ALSUTRECHT.epochRejections.InterpReport    = report;
EEG.ALSUTRECHT.epochRejections.interpEpochs    = length(report.listFixed);

num_trials(3) = EEG.trials;
fprintf('Retained %d/%d trials after Stage B interpolation/rejection.\n', num_trials(3), num_trials(1));

% =============================================
% Note the final number of trials
% =============================================
num_trials(4) = EEG.trials;

% =========================================================================
% Log
EEG.ALSUTRECHT.epochRejections.epoch_initial     = num_trials(1);
EEG.ALSUTRECHT.epochRejections.epoch_stage_1     = num_trials(2);
EEG.ALSUTRECHT.epochRejections.epoch_stage_2     = num_trials(3);
EEG.ALSUTRECHT.epochRejections.epoch_final       = num_trials(4);
EEG.ALSUTRECHT.epochRejections.mask_eeglab       = mask_eeglab;
EEG.ALSUTRECHT.epochRejections.mask_amplitude    = mask_amplitude;
EEG.ALSUTRECHT.epochRejections.mask_muscle_env   = mask_muscle_env;
EEG.ALSUTRECHT.epochRejections.mask_muscle_slope = mask_muscle_slope;
EEG.ALSUTRECHT.epochRejections.reject_stats      = reject_stats;
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


function mask_amplitude = detect_high_amplitude(EEG, cfg)
freq_stop = [2 40]; % Hz
freq_sample = EEG.srate;

% Prototype order 2 -> 4th-order IIR -> 8th-order effective in filtfilt
% Provides strong stopband attenuation without long settling times
[b, a] = butter(2, freq_stop / (freq_sample / 2), 'stop');
assert(isstable(b, a), 'Bandstop filter unstable.');

% Extract channels [channels x timepoints x epochs]
chan_mask = strcmp({EEG.chanlocs.type}, 'EEG');
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
edge_trim_pts = max(1, round(0.050 * freq_sample));
eval_idx = (1 + edge_trim_pts) : (n_pnts - edge_trim_pts);

% 4. Evaluate maximum absolute amplitude across valid time points and channels
% Shape: [channels x epochs]
max_per_chan_trial = squeeze(max(abs(data_filt(eval_idx, :, :)), [], 1));

% 5. Trial rejection mask: flag trial if any channel exceeds threshold
mask_amplitude = any(max_per_chan_trial > cfg.epoch.amplitude_max, 1)';

fprintf('%d/%d trials marked for rejection.\n', sum(mask_amplitude), n_trials);

end

function mask_envelope = detect_high_muscle(EEG, cfg)
freq_muscle = 70; % Hz
freq_sample = EEG.srate;

% Check the lowpass fitlering settings
cfg_filt = extract_lpfilt(EEG.ALSUTRECHT.subject.task, cfg);
if isempty(cfg_filt.lp)
    cfg_filt.lp(1) = Inf;
end

% Only run EMG envelope detection if data bandwidth permits frequencies >70 Hz
if (cfg_filt.lp(1) - freq_muscle) > 20
    chan_mask = strcmp({EEG.chanlocs.type}, 'EEG');
    data_tmp = double(EEG.data(chan_mask, :, :));

    % Concatenate channels across continuous time
    [n_ch, n_pnts, n_trials] = size(data_tmp);
    data_tmp = reshape(data_tmp, n_ch, []);

    mask_envelope = detect_emg_envelopes(data_tmp, freq_sample, freq_muscle);
    mask_envelope = reshape(mask_envelope, [n_pnts, n_trials]);
    mask_envelope = any(mask_envelope, 1)';

    fprintf('%d/%d trials marked for rejection by EMG envelope.\n', sum(mask_envelope), n_trials);
else
    fprintf('Bypassing EMG envelope check (Data low-pass filtered at %d Hz; requires >%d Hz).\n', cfg_filt.lp(1), freq_muscle + 20);
    mask_envelope = false(size(EEG.data, 3), 1);
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
end

function [diag_stats, fh] = check_rejected_trials(EEG, rej_idx, cfg)
% CHECK_REJECTED_TRIALS
% Audits trials flagged for rejection by EEGLAB statistical tools
% (e.g. pop_jointprob, pop_rejkurt) to determine whether pruning targets
% genuine artefacts (EMG/drift) or physiological rhythms (alpha/beta bursts).
%
% Renders a 2x2 diagnostic figure:
%   1. Rejection Timeline (flagged distribution, rolling rate, EO/EC split)
%   2. Spectral Density Overlay (Kept vs Flagged, 1-50 Hz)
%   3. Surplus RMS Broadband Topography (Flagged - Kept)
%   4. Surplus RMS High-Frequency Topography (>30 Hz)

% -------------------------------------------------------------------------
% 1. Parameter Validation & Channel Extraction
% -------------------------------------------------------------------------
if nargin < 3, cfg = struct(); end
if ~isfield(cfg, 'hf_cutoff'), cfg.hf_cutoff = 30;   end % Hz
if ~isfield(cfg, 'f_max'),     cfg.f_max     = 50;   end % Hz
if ~isfield(cfg, 'do_plot'),   cfg.do_plot   = true; end
if ~isfield(cfg, 'sub_id'),    cfg.sub_id    = '';   end

ch_idx = find(strcmpi({EEG.chanlocs.type}, 'EEG'));
if isempty(ch_idx), ch_idx = 1:EEG.nbchan; end
n_chans = length(ch_idx);
n_trials = EEG.trials;
fs = EEG.srate;

% -------------------------------------------------------------------------
% 2. Resolve Rejected Trial Indices & Task Transition (EO -> EC)
% -------------------------------------------------------------------------
if nargin < 2 || isempty(rej_idx)
    rej_mask = false(1, n_trials);
    if isfield(EEG, 'reject') && ~isempty(EEG.reject)
        fnames = fieldnames(EEG.reject);
        for f = 1:numel(fnames)
            val = EEG.reject.(fnames{f});
            if islogical(val) && length(val) == n_trials
                rej_mask = rej_mask | val;
            end
        end
    end
else
    if islogical(rej_idx)
        rej_mask = rej_idx(:)';
    else
        rej_mask = false(1, n_trials);
        rej_mask(rej_idx) = true;
    end
end

n_rej  = sum(rej_mask);
n_kept = sum(~rej_mask);

% Default: split at midpoint
split_idx   = floor(n_trials / 2) + 1;
label_h1    = '1st Half';
label_h2    = '2nd Half';
split_label = 'Midpoint';

% Detect RS transition from Eyes Open (EO) to Eyes Closed (EC)
is_rs = isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'subject') && ...
    isfield(EEG.ALSUTRECHT.subject, 'task') && strcmpi(EEG.ALSUTRECHT.subject.task, 'RS');

if is_rs && isfield(EEG, 'event') && ~isempty(EEG.event)
    ev_types = cellfun(@num2str, {EEG.event.type}, 'UniformOutput', false);

    % Find EC1 or the first event starting with EC
    ec_ev_idx = find(strcmpi(ev_types, 'EC1'), 1);
    if isempty(ec_ev_idx)
        ec_ev_idx = find(startsWith(ev_types, 'EC', 'IgnoreCase', true), 1);
    end

    if ~isempty(ec_ev_idx)
        % In epoched datasets, .epoch resolves the exact trial regardless of event counts
        if isfield(EEG.event, 'epoch') && ~isempty(EEG.event(ec_ev_idx).epoch)
            split_idx = EEG.event(ec_ev_idx).epoch;
        else
            % Fallback mapping via event-to-trial ratio
            ev_ratio  = length(EEG.event) / n_trials;
            split_idx = max(1, min(n_trials, ceil(ec_ev_idx / ev_ratio)));
        end
        label_h1    = 'Eyes Open (EO)';
        label_h2    = 'Eyes Closed (EC)';
        split_label = 'EO -> EC Transition';
    end
end

% Compute split statistics
idx_h1 = 1:(split_idx - 1);
idx_h2 = split_idx:n_trials;

rej_h1 = sum(rej_mask(idx_h1));
rej_h2 = sum(rej_mask(idx_h2));
pct_h1 = (rej_h1 / max(1, length(idx_h1))) * 100;
pct_h2 = (rej_h2 / max(1, length(idx_h2))) * 100;

fprintf('--- Rejection Audit: %d Flagged / %d Total (%.1f%%) ---\n', ...
    n_rej, n_trials, (n_rej / n_trials) * 100);
fprintf('    - %s (Trials 1-%d): %d flagged (%.1f%%)\n', ...
    label_h1, split_idx - 1, rej_h1, pct_h1);
fprintf('    - %s (Trials %d-%d): %d flagged (%.1f%%)\n', ...
    label_h2, split_idx, n_trials, rej_h2, pct_h2);

if n_rej == 0 || n_kept == 0
    warning('Cannot run comparison: %d flagged, %d retained.', n_rej, n_kept);
    diag_stats = struct();
    fh = [];
    return;
end

kept_mask = ~rej_mask;

% -------------------------------------------------------------------------
% 3. Diagnostic 1: Spectral Profile (Welch PSD)
% -------------------------------------------------------------------------
n_pnts = size(EEG.data, 2);
n_fft  = min(n_pnts, round(2.0 * fs));
win    = hann(n_fft);
n_ovlp = round(n_fft / 2);

data_all = double(EEG.data(ch_idx, :, :));
data_all = data_all - mean(data_all, 2);

[~, f_axis] = pwelch(data_all(1, :, 1), win, n_ovlp, n_fft, fs);
n_freqs = length(f_axis);

psd_matrix = zeros(n_chans, n_freqs, n_trials);
for tr = 1:n_trials
    psd_matrix(:, :, tr) = pwelch(data_all(:, :, tr)', win, n_ovlp, n_fft, fs)';
end

psd_kept_mean = mean(mean(psd_matrix(:, :, kept_mask), 3), 1);
psd_rej_mean  = mean(mean(psd_matrix(:, :, rej_mask), 3), 1);

idx_alpha = f_axis >= 8  & f_axis <= 13;
idx_hf    = f_axis >= 35 & f_axis <= min(cfg.f_max, 45);

alpha_power_rej = trapz(f_axis(idx_alpha), psd_rej_mean(idx_alpha));
hf_power_rej    = trapz(f_axis(idx_hf), psd_rej_mean(idx_hf));
ratio_alpha_hf  = alpha_power_rej / max(hf_power_rej, eps);

% -------------------------------------------------------------------------
% 4. Diagnostic 2: Spatial Topographies of Surplus RMS
% -------------------------------------------------------------------------
rms_per_trial = squeeze(sqrt(mean(data_all.^2, 2)));
rms_kept_ch   = mean(rms_per_trial(:, kept_mask), 2);
rms_rej_ch    = mean(rms_per_trial(:, rej_mask), 2);
rms_diff      = rms_rej_ch - rms_kept_ch;

[b_hf, a_hf] = butter(4, cfg.hf_cutoff / (fs / 2), 'high');
data_flat    = reshape(permute(data_all, [2, 1, 3]), n_pnts, []);
data_hf_flat = filtfilt(b_hf, a_hf, data_flat);
data_hf      = permute(reshape(data_hf_flat, n_pnts, n_chans, n_trials), [2, 1, 3]);

rms_hf_per_trial = squeeze(sqrt(mean(data_hf.^2, 2)));
rms_hf_kept_ch   = mean(rms_hf_per_trial(:, kept_mask), 2);
rms_hf_rej_ch    = mean(rms_hf_per_trial(:, rej_mask), 2);
rms_hf_diff      = rms_hf_rej_ch - rms_hf_kept_ch;

% -------------------------------------------------------------------------
% 5. Pack Output Statistics
% -------------------------------------------------------------------------
diag_stats.f_axis         = f_axis;
diag_stats.psd_kept_mean  = psd_kept_mean;
diag_stats.psd_rej_mean   = psd_rej_mean;
diag_stats.rms_diff       = rms_diff;
diag_stats.rms_hf_diff    = rms_hf_diff;
diag_stats.ratio_alpha_hf = ratio_alpha_hf;
diag_stats.n_rej          = n_rej;
diag_stats.n_kept         = n_kept;
diag_stats.rej_mask       = rej_mask;
diag_stats.split_idx      = split_idx;
diag_stats.label_h1       = label_h1;
diag_stats.label_h2       = label_h2;
diag_stats.pct_h1         = pct_h1;
diag_stats.pct_h2         = pct_h2;
diag_stats.n_rej_h1       = rej_h1;
diag_stats.n_rej_h2       = rej_h2;

% -------------------------------------------------------------------------
% 6. Visualisation (2 Rows x 2 Columns)
% -------------------------------------------------------------------------
fh = [];
if cfg.do_plot
    fh = figure('Color', 'w', 'Position', [100, 100, 1200, 750]);

    % Panel 1: Rejection Timeline & Condition Transition
    subplot(2, 2, 1);
    hold on;
    xline(split_idx - 0.5, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4, ...
        'DisplayName', split_label);

    idx_flagged = find(rej_mask);
    if ~isempty(idx_flagged)
        stem(idx_flagged, ones(size(idx_flagged)), 'Color', [0.85 0.20 0.20], ...
            'MarkerFaceColor', [0.85 0.20 0.20], 'LineWidth', 1.1, 'MarkerSize', 4, ...
            'DisplayName', sprintf('Flagged (N = %d)', n_rej));
    end

    % Overlaid rolling rejection rate (window ~ 10% of total trials)
    win_mov  = max(5, round(n_trials * 0.10));
    rate_mov = movmean(double(rej_mask), win_mov);
    plot(1:n_trials, rate_mov, 'Color', [0.15 0.50 0.80], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Rolling Rate (win = %d)', win_mov));

    grid on; box off;
    xlim([1, n_trials]);
    ylim([0, 1.15]);
    set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1], 'YTickLabel', {'0%', '25%', '50%', '75%', '100%'});
    xlabel('Trial / Epoch Index', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Flagged Trials & Rate', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Rejection Timeline (%s: %.1f%% | %s: %.1f%%)', ...
        label_h1, pct_h1, label_h2, pct_h2), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'Box', 'off', 'FontSize', 8);

    % Panel 2: Spectral Density Overlay
    subplot(2, 2, 2);
    idx_plot = f_axis >= 1 & f_axis <= cfg.f_max;
    plot(f_axis(idx_plot), 10*log10(psd_kept_mean(idx_plot)), 'Color', [0.15 0.50 0.80], ...
        'LineWidth', 2.0, 'DisplayName', sprintf('Kept (N = %d)', n_kept)); hold on;
    plot(f_axis(idx_plot), 10*log10(psd_rej_mean(idx_plot)), 'Color', [0.85 0.20 0.20], ...
        'LineWidth', 2.0, 'DisplayName', sprintf('Flagged (N = %d)', n_rej));
    grid on; box off;
    xlim([1, cfg.f_max]);
    xlabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Power (10\cdotlog_{10} \muV^2/Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('PSD Profile (Alpha/HF Ratio: %.2f)', ratio_alpha_hf), ...
        'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'Box', 'off');

    % Symmetric colour limits
    max_d_tot = max(abs(rms_diff));
    if max_d_tot == 0, max_d_tot = 1; end
    clim_tot = [-ceil(max_d_tot * 1.1), ceil(max_d_tot * 1.1)];

    max_d_hf = max(abs(rms_hf_diff));
    if max_d_hf == 0, max_d_hf = 1; end
    clim_hf = [-ceil(max_d_hf * 1.1), ceil(max_d_hf * 1.1)];

    % Panel 3: Broadband Surplus Topography
    ax3 = subplot(2, 2, 3);
    mytopoplot(rms_diff, [], 'Surplus RMS: Broadband (Flagged - Kept)', ax3, clim_tot);
    try
        colormap(ax3, brewermap([], '*RdBu'));
    catch
        colormap(ax3, jet);
    end
    clim(ax3, clim_tot);
    cb3 = colorbar(ax3);
    ylabel(cb3, '\Delta \muV', 'FontSize', 9, 'FontWeight', 'bold');

    % Panel 4: High-Frequency Surplus Topography (>30 Hz)
    ax4 = subplot(2, 2, 4);
    mytopoplot(rms_hf_diff, [], sprintf('Surplus RMS: >%d Hz (Flagged - Kept)', cfg.hf_cutoff), ax4, clim_hf);
    try
        colormap(ax4, brewermap([], '*RdBu'));
    catch
        colormap(ax4, jet);
    end
    clim(ax4, clim_hf);
    cb4 = colorbar(ax4);
    ylabel(cb4, '\Delta \muV', 'FontSize', 9, 'FontWeight', 'bold');

    if ~isempty(cfg.sub_id)
        sgtitle(sprintf('Rejection Diagnostic Audit: %s', cfg.sub_id), ...
            'FontSize', 13, 'FontWeight', 'bold');
    end
end

end