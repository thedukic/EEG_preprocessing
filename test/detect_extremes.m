function [EEG, flag_exclude] = detect_extremes(EEG, cfg)
% DETECT_EXTREMES Screen and reject extreme EMG and postural artifacts across blocks
%
% Usage:
%   [EEG, bad_masks] = detect_extremes(EEG, cfg);
%
% Inputs:
%   EEG - EEGLAB dataset struct or array of structs (e.g. EEG(1:6))
%   cfg - Configuration structure with optional fields:
%         .z_emg_high      : Robust Z-score seed trigger for 50-100 Hz EMG (default: 5.0)
%         .z_emg_low       : Robust Z-score boundary floor for hysteresis (default: 3.0)
%         .min_emg_uV_high : Absolute physical floor in uV for EMG seed (default: 10.0)
%         .min_emg_uV_low  : Absolute physical floor in uV for boundary (default: 6.0)
%         .min_dur_sec     : Minimum continuous duration of raw EMG burst before dilation (default: 0.05)
%         .z_swing         : Robust Z-score threshold for dominant covariance swings (default: 6.5)
%         .win_sec         : Sliding window length in seconds for covariance estimation (default: 1.0)
%         .min_mad_cov     : Floor for covariance MAD to protect blinks in ultra-clean data (default: 0.35)
%         .pad_sec         : Dilation padding applied symmetrically on each side (default: 0.10)
%         .min_good_sec    : Minimum duration in seconds of clean data to retain (default: 0.0)
%                            Clean segments <= this duration are absorbed and padded as bad.
%         .hp_filter       : Apply zero-phase 0.3 Hz highpass filter for baseline stability (default: true)
%         .plot            : Open interactive eegplot per block (default: true)
%         .opt_figure      : Visibility of summary figure ('on' or 'off', default: 'on')

% -------------------------------------------------------------
% 0. INITIALISE ALL DEFAULTS AT TOP LEVEL
% -------------------------------------------------------------
if nargin < 2, cfg = struct(); end

% Covariance and swing detection (Gross mechanical/cable pulls only)
if ~isfield(cfg, 'z_swing'),      cfg.z_swing      = 6.5;  end
if ~isfield(cfg, 'win_sec'),      cfg.win_sec      = 1.0;  end
if ~isfield(cfg, 'min_mad_cov'),  cfg.min_mad_cov  = 0.35; end
if ~isfield(cfg, 'pad_sec'),      cfg.pad_sec      = 0.10; end

% EMG Thresholds: STRICTLY EXTREME BURSTS
if ~isfield(cfg, 'z_emg_high'),      cfg.z_emg_high      = 5.0;  end
if ~isfield(cfg, 'min_emg_uV_high'), cfg.min_emg_uV_high = 12.0; end % Ignore minor twitches

if ~isfield(cfg, 'z_emg_low'),       cfg.z_emg_low       = 3.0;  end
if ~isfield(cfg, 'min_emg_uV_low'),  cfg.min_emg_uV_low  = 6.0;  end % Stop greedy expansion!
if ~isfield(cfg, 'min_dur_sec'),     cfg.min_dur_sec     = 0.05; end

% STOP the script from eating clean data between bursts
if ~isfield(cfg, 'min_good_sec'),    cfg.min_good_sec    = 2.00; end % Do not absorb clean gaps
if ~isfield(cfg, 'pad_sec'),         cfg.pad_sec         = 0.10; end % Keep padding minimal

% Preprocessing and interface options
if ~isfield(cfg, 'hp_filter'),   cfg.hp_filter   = true;  end
if ~isfield(cfg, 'plot'),        cfg.plot        = false; end
if ~isfield(cfg, 'opt_figure'),  cfg.opt_figure  = 'on';  end

num_block    = numel(EEG);
mask_bad     = cell(1, num_block);
interval_bad = cell(1, num_block);

% -------------------------------------------------------------
% 1. DETECTION LOOP ACROSS BLOCKS (PRE-REJECTION)
% -------------------------------------------------------------
perc_segment = NaN(num_block, 1);
for i_block = 1:num_block
    block_name = '';
    if isfield(EEG(i_block), 'setname') && ~isempty(EEG(i_block).setname)
        block_name = sprintf(' [%s]', EEG(i_block).setname);
    end
    fprintf('\n--- Screening Block %d of %d%s ---\n', i_block, num_block, block_name);

    % Call helper detection function on original un-cut block
    [mask_bad{i_block}, interval_bad{i_block}, plot_data] = detect_extremes_single(EEG(i_block), cfg);

    % Save mask and intervals into dataset hierarchy
    EEG(i_block).ALSUTRECHT.bad_segments.mask_bad     = mask_bad;
    EEG(i_block).ALSUTRECHT.bad_segments.interval_bad = interval_bad;

    num_segment = size(interval_bad{i_block}, 1);
    perc_segment(i_block) = sum(mask_bad{i_block}) / EEG(i_block).pnts;
    fprintf('Detected %d bad segments (%.2f%% of data).\n', num_segment, perc_segment(i_block) * 100);

    % Optional per-block interactive scrolling viewer
    if cfg.plot && num_segment > 0
        plot_mask(EEG(i_block), plot_data, mask_bad{i_block});
    end
end

% Weighted mean across all session samples:
perc_bad = sum(cellfun(@sum, mask_bad)) / sum([EEG.pnts]);
flag_exclude = perc_bad > 0.50;

% -------------------------------------------------------------
% 2. MULTI-BLOCK SUMMARY VISUALISATION
% -------------------------------------------------------------
report_badsegments(EEG, interval_bad, 'bad_segments', cfg.opt_figure);

% -------------------------------------------------------------
% 3. REMOVAL LOOP ACROSS BLOCKS
% -------------------------------------------------------------
for i_block = 1:num_block
    intervals = interval_bad{i_block};
    if ~isempty(intervals)
        fprintf('Rejecting %d segment(s) in Block %d...\n', size(intervals, 1), i_block);
        EEG(i_block) = eeg_eegrej(EEG(i_block), intervals);
        EEG(i_block) = eeg_checkset(EEG(i_block));
    end

    EEG(i_block).ALSUTRECHT.bad_segments.perc_bad = perc_bad;
end

% if num_block == 1
%     bad_masks = bad_masks{1};
% end

end


% =============================================================
% HELPER: DETECTION ON A SINGLE BLOCK
% =============================================================
function [bad_mask, intervals, plot_data] = detect_extremes_single(EEG, cfg)

srate     = EEG.srate;
n_samples = EEG.pnts;

% 1. Optional 0.3 Hz High-pass Filter for baseline stability
if cfg.hp_filter
    [b_hp, a_hp] = butter(2, 0.3 / (srate / 2), 'high');
    plot_data = filtfilt(b_hp, a_hp, double(EEG.data)')';
else
    plot_data = double(EEG.data);
end

% 2. Isolate EEG Channels Only
if isfield(EEG.chanlocs, 'type') && any(~cellfun(@isempty, {EEG.chanlocs.type}))
    channel_mask = strcmpi({EEG.chanlocs.type}, 'EEG');
    if ~any(channel_mask)
        warning('No channels found with type == "EEG". Using all channels.');
        channel_mask = true(1, EEG.nbchan);
    end
else
    channel_mask = true(1, EEG.nbchan);
end
eeg_data = plot_data(channel_mask, :);

% -------------------------------------------------------------
% 1. FOCAL HIGH-FREQUENCY EMG DETECTOR (Two-Pass + Hysteresis)
% -------------------------------------------------------------
[b_hf, a_hf]   = butter(4, [50 100] / (srate / 2), 'bandpass');
[b_env, a_env] = butter(2, 2 / (srate / 2), 'low');

hf_data = filtfilt(b_hf, a_hf, eeg_data')';
env_sq  = filtfilt(b_env, a_env, (hf_data.^2)')';
hf_env  = sqrt(max(env_sq, 0));

% 95th percentile channel isolates focal temporal/cranial spikes
focal_emg_trace = prctile(hf_env, 95, 1);

% PASS 1A: Preliminary baseline to identify uncontaminated samples
init_med = median(focal_emg_trace);
init_mad = 1.4826 * mad(focal_emg_trace, 1);
clean_samples = focal_emg_trace < (init_med + 3.0 * init_mad);

if sum(clean_samples) < (0.2 * n_samples)
    clean_samples = true(1, n_samples);
end

% PASS 1B: True baseline stats
med_emg = median(focal_emg_trace(clean_samples));
mad_emg = 1.4826 * mad(focal_emg_trace(clean_samples), 1);
z_emg   = (focal_emg_trace - med_emg) / (mad_emg + eps);

% PASS 1C: Hysteresis Thresholding (Seed & Grow)
seed_raw = (z_emg > cfg.z_emg_high) & (focal_emg_trace > cfg.min_emg_uV_high);
boundary_mask = (z_emg > cfg.z_emg_low) & (focal_emg_trace > cfg.min_emg_uV_low);

% Diagnostic Print
active_trigger_uV = max(cfg.min_emg_uV_high, med_emg + cfg.z_emg_high * mad_emg);
fprintf('   [Diag] Baseline: %.2f uV | True MAD: %.2f | Block Max: %.2f uV | ACTIVE TRIGGER: %.2f uV\n', ...
    med_emg, mad_emg, max(focal_emg_trace), active_trigger_uV);

% Eliminate single-sample digital spikes
min_samples = round(cfg.min_dur_sec * srate);
seed_mask   = bwareaopen(seed_raw, min_samples);

if any(seed_mask)
    emg_mask = imreconstruct(seed_mask, boundary_mask);
else
    emg_mask = false(1, n_samples);
end

% Bridge brief motor unit firing flutter gaps (< 150 ms)
emg_mask = imclose(emg_mask, strel('line', round(0.15 * srate), 0));

% -------------------------------------------------------------
% 2. GROSS POSTURAL & CABLE SWING DETECTOR
% -------------------------------------------------------------
win_len      = round(cfg.win_sec * srate);
step_size    = round(win_len / 4);
win_starts   = 1:step_size:(n_samples - win_len);
n_wins       = length(win_starts);
log_lambda_1 = zeros(1, n_wins);

for w = 1:n_wins
    idx = win_starts(w) : (win_starts(w) + win_len - 1);
    X = eeg_data(:, idx);
    X = X - mean(X, 2);

    C = (X * X') / win_len;
    log_lambda_1(w) = log(real(eigs(C, 1, 'largestabs')) + eps);
end

med_lambda = median(log_lambda_1);
mad_lambda = 1.4826 * mad(log_lambda_1, 1);
mad_lambda = max(mad_lambda, cfg.min_mad_cov);

z_swing = (log_lambda_1 - med_lambda) / mad_lambda;

bad_swing_wins = z_swing > cfg.z_swing;
swing_mask     = false(1, n_samples);
for w = 1:n_wins
    if bad_swing_wins(w)
        swing_mask(win_starts(w) : win_starts(w) + win_len - 1) = true;
    end
end

swing_mask = imclose(swing_mask, strel('line', round(0.5 * srate), 0));

% -------------------------------------------------------------
% 3. MERGE MASKS, DILATE & ABSORB SHORT GOOD SEGMENTS
% -------------------------------------------------------------
raw_bad_mask = emg_mask | swing_mask;

if any(raw_bad_mask)
    closed_mask = imclose(raw_bad_mask, strel('line', round(0.3 * srate), 0));
    pad_samples = round(cfg.pad_sec * srate);
    bad_mask    = imdilate(closed_mask, strel('line', pad_samples * 2, 0));
else
    bad_mask = false(1, n_samples);
end

% Safely bypass fragmentation absorption if min_good_sec is 0
if cfg.min_good_sec > 0
    min_good_samples = round(cfg.min_good_sec * srate);
    good_mask_valid  = bwareaopen(~bad_mask, min_good_samples + 1);
    bad_mask         = ~good_mask_valid;
end

% -------------------------------------------------------------
% 4. EXTRACT FINAL BOUNDARIES
% -------------------------------------------------------------
diff_mask = diff([0, bad_mask, 0]);
starts    = find(diff_mask == 1);
stops     = find(diff_mask == -1) - 1;

intervals = [starts(:), stops(:)];

end


% =============================================================
% HELPER: INTERACTIVE SCROLLING VIEWER
% =============================================================
function plot_mask(EEG, plot_data, burst_mask)
diff_mask = diff([0, burst_mask(:)', 0]);
onsets    = find(diff_mask == 1);
offsets   = find(diff_mask == -1) - 1;

n_all_chans = size(plot_data, 1);
n_bursts    = numel(onsets);

winrej = zeros(n_bursts, 5 + n_all_chans);
winrej(:, 1)     = onsets(:);
winrej(:, 2)     = offsets(:);
winrej(:, 3:5)   = repmat([0.9 0.4 0.4], n_bursts, 1);
winrej(:, 6:end) = 1;

fig_title = sprintf('Extreme Artifacts (Count: %d)', n_bursts);
if isfield(EEG, 'setname') && ~isempty(EEG.setname)
    fig_title = sprintf('%s - %s', EEG.setname, fig_title);
end

eegplot(plot_data, 'srate', EEG.srate, ...
    'title', fig_title, ...
    'eloc_file', EEG.chanlocs, ...
    'winrej', winrej, ...
    'winlength', 10, ...
    'spacing', 100, ...
    'color', 'off');
end