function check_power_cleaning(EEG_old, EEG_new, ics_to_remove, figure_tag, cfg)
% CHECK_SPECTRAL_POWER compares PSD before and after cleaning.
% Automatically selects the Region of Interest (ROI) and frequency band
% based on the task defined in EEG_old.ALSUTRECHT.subject.task.

% 1. Determine parameters based on task type
if isfield(EEG_old.ALSUTRECHT.subject, 'task')
    task = EEG_old.ALSUTRECHT.subject.task;
else
    warning('Task not defined in EEG structure. Defaulting to Global.');
    task = 'UNKNOWN';
end

switch upper(task)
    case 'RS'
        target_labels = cfg.roi.rs_alpha;
        roi_name      = 'Occipital';
        target_band   = [7 13];
        band_name     = 'Alpha Band (7-13 Hz)';

    case 'MMN'
        target_labels = cfg.roi.mmn;
        roi_name      = 'Frontocentral-MMN';
        target_band   = [4 8];
        band_name     = 'Theta Band (4-8 Hz)';

    case 'SART'
        target_labels = cfg.roi.sart_p300;
        roi_name      = 'Centroparietal-P300';
        target_band   = [4 8];
        band_name     = 'Theta/P300 Band (4-8 Hz)';

    case 'MT'
        target_labels = cfg.roi.motor_all;
        roi_name      = 'Sensorimotor';
        target_band   = [10 30];
        band_name     = 'Mu/Beta Band (10-30 Hz)';

    otherwise
        target_labels = {};
        roi_name      = 'Global';
        target_band   = [];
        band_name     = '';
end

% 2. Identify target channels
all_chan_labels = {EEG_old.chanlocs.labels};
if isempty(target_labels)
    chan_idx = 1:length(all_chan_labels);
else
    chan_idx = find(ismember(all_chan_labels, target_labels));
    if isempty(chan_idx)
        warning('Specified target channels not found. Defaulting to all channels.');
        chan_idx = 1:length(all_chan_labels);
        roi_name = 'Global';
    end
end

% 3. Extract raw data and back-project cleaned data
if isempty(EEG_new) && ~isempty(ics_to_remove)
    assert(isfield(EEG_old.ALSUTRECHT, 'ica'));
    is_ica = true;

    if ~isempty(EEG_old.ALSUTRECHT.ica.icachansind)
        ic_ch_idx = EEG_old.ALSUTRECHT.ica.icachansind;
    else
        ic_ch_idx = 1:size(EEG_old.data, 1);
    end
elseif ~isempty(EEG_new) && isempty(ics_to_remove)
    is_ica = false;
    ic_ch_idx = 1:size(EEG_old.data, 1);
else
    error('Invalid input combination: Provide either EEG_new OR ics_to_remove.');
end

raw_all_2d = reshape(EEG_old.data(ic_ch_idx, :, :), length(ic_ch_idx), []);

% Extract the clean data
if ~is_ica
    % Simply use the cleaned data that is passed on
    clean_all_2d = EEG_new.data(ic_ch_idx, :);

elseif is_ica
    % Calculate component activations and zero out removed ICs
    W = EEG_old.ALSUTRECHT.ica.icaweights * EEG_old.ALSUTRECHT.ica.icasphere;
    activations = W * raw_all_2d;
    activations(ics_to_remove, :) = 0;

    % Reconstruct channel space
    Winv = EEG_old.ALSUTRECHT.ica.icawinv;
    clean_all_2d = Winv * activations;
end

% Isolate target channels for the ROI
[~, loc_in_ica] = ismember(chan_idx, ic_ch_idx);
loc_in_ica = loc_in_ica(loc_in_ica > 0);

data_raw_roi   = raw_all_2d(loc_in_ica, :);
data_clean_roi = clean_all_2d(loc_in_ica, :);

% 4. Compute Power Spectral Density (PSD) using Welch's method
fs = EEG_old.srate;
win_len = 2 * fs; % 2-second Hamming window
noverlap = win_len / 2;

% Average across the selected ROI channels
sig_raw   = mean(data_raw_roi, 1);
sig_clean = mean(data_clean_roi, 1);

[psd_raw, f]   = pwelch(sig_raw, win_len, noverlap, [], fs);
[psd_clean, ~] = pwelch(sig_clean, win_len, noverlap, [], fs);

% Convert to dB power
psd_raw_db   = 10 * log10(psd_raw);
psd_clean_db = 10 * log10(psd_clean);

% 5. Plot comparison
fh = figure('Color', 'w', 'Position', [100, 100, 750, 480], 'Visible', cfg.figure.visible);
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(t);

% Plot raw vs cleaned spectra
plot(f, psd_raw_db, 'Color', [0.8, 0.25, 0.25], 'LineWidth', 1.5, 'DisplayName', 'Pre-Cleaning');
hold on;
plot(f, psd_clean_db, 'Color', [0.1, 0.6, 0.3], 'LineWidth', 1.8, 'DisplayName', 'Post-Cleaning');

% Highlight Specific Frequency Band (if requested)
f_sub = (f >= 1 & f <= 45);

min_1 = min(psd_raw_db(f_sub));
min_2 = min(psd_clean_db(f_sub));
max_1 = max(psd_raw_db(f_sub));
max_2 = max(psd_clean_db(f_sub));

y_lims = [min(min_1, min_2) - 2, max(max_1, max_2) + 3];

if ~isempty(target_band)
    patch([target_band(1) target_band(2) target_band(2) target_band(1)], ...
        [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', band_name);

    % Restack lines over patch
    uistack(findobj(ax, 'Type', 'line'), 'top');
end

% Styling
xlim([1, 45]);
ylim(y_lims);
grid on;
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power Spectral Density (dB/Hz)', 'FontSize', 12, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Dynamic Title Calculation
% -------------------------------------------------------------------------
if ~isempty(target_band)
    idx_band = (f >= target_band(1) & f <= target_band(2));
    diff_label = 'Target Band';
else
    idx_band = (f >= 1 & f <= 45);
    diff_label = 'Broadband';
end

% 1. Calculate dB difference (Logarithmic)
power_diff_db = mean(psd_clean_db(idx_band) - psd_raw_db(idx_band));

% 2. Calculate percentage difference (Linear)
% Use the linear PSD values, average them across the band, then find the % change
mean_p_raw   = mean(psd_raw(idx_band));
mean_p_clean = mean(psd_clean(idx_band));
power_diff_pct = ((mean_p_clean - mean_p_raw) / mean_p_raw) * 100;

% 3. Format the title string
% if is_ica
%     num_removed = sum(ics_to_remove);
%     title_str = sprintf('%s PSD Check (N = %d ICs Removed) | %s \\Delta: %.2f dB (%.1f%%)', ...
%         roi_name, num_removed, diff_label, power_diff_db, power_diff_pct);
% else
%     title_str = sprintf('%s PSD Check | %s \\Delta: %.2f dB (%.1f%%)', ...
%         roi_name, diff_label, power_diff_db, power_diff_pct);
% end
title_str = sprintf('%s PSD Check | %s \\Delta: %.2f dB (%.1f%%)', ...
    roi_name, diff_label, power_diff_db, power_diff_pct);

title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
% -------------------------------------------------------------------------

legend('Location', 'northeast', 'Box', 'on');
hold off;

% Save
save_figure(fh, EEG_old.ALSUTRECHT.subject.figures, [EEG_old.ALSUTRECHT.subject.id '_power_post_' figure_tag], [20 20/1.6]);

end