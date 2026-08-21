function EEG = estimate_electrodeoffsets(EEG, cfg)
% ESTIMATE_ELECTRODE_OFFSETS calculates the DC offset for BioSemi electrodes.
%
% Inputs:
%   EEG          : Vector of EEGLAB structures (1 x num_blocks).
%                  Data MUST be raw. No high-pass filtering, detrending,
%                  or average referencing can be applied prior to this step.
%                  Assumes input values are in microvolts (uV).
%   threshold_mV : numeric scalar
%                  The threshold in millivolts to flag poor contact.
%                  BioSemi typically recommends +/- 40 mV.
%
% Outputs:
%   offsets_mV   : numeric matrix (channels x blocks) of estimated offsets.
%   bad_channels : cell array (1 x blocks) of channel indices exceeding the threshold.

% 50 is max, but 25 is ideal
threshold_mV = 50;

fprintf('\n================================\n');
fprintf('Estimating electrode offsets\n');
fprintf('================================\n');

num_block = length(EEG);
num_channels = EEG(1).nbchan;

% Preallocate outputs for performance
offsets_mV = NaN(num_channels, num_block);
bad_electrodes = cell(1, num_block);

for i_block = 1:num_block
    % 1. Extract the DC component by averaging over time
    % Access the data for the specific block using EEG(i_block)
    offsets_uV = mean(EEG(i_block).data, 2);

    % 2. Convert from microvolts to millivolts
    % Assign the column vector to the specific block column in the matrix
    offsets_mV(:, i_block) = offsets_uV / 1000;

    % 3. Identify channels that exceed the absolute threshold
    % Store in a cell array because the number of indices can vary per block
    bad_electrodes{i_block} = find(abs(offsets_mV(:, i_block)) > threshold_mV);
end

% Display summary
fprintf('Calculated offsets for %d channels across %d blocks.\n', num_channels, num_block);
fprintf('Average absolute electrode offset: %1.1f mV\n', mean(abs(offsets_mV), "all"));

% Find unique bad channels across all blocks for the summary warning
bad_electrodes_all_indx = unique(cat(1, bad_electrodes{:}));
bad_electrodes_all_indx = bad_electrodes_all_indx(bad_electrodes_all_indx <= 128);

channel_mask = strcmp({EEG(1).chanlocs.type}, 'EEG');
channel_labels = {EEG(1).chanlocs(channel_mask).labels};
bad_electrodes_all_label = channel_labels(bad_electrodes_all_indx);

if ~isempty(bad_electrodes_all_label)
    fprintf('Unusually high offsets (> %d mV; unused/broken electrodes?) found (N = %d):\n', threshold_mV, length(bad_electrodes_all_label));
    fprintf('--------------------------------------------------\n');

    % Extract the specific labels and values for the strange channels
    max_offset = max(abs(offsets_mV), [], 2);
    strange_offsets = max_offset(bad_electrodes_all_indx);

    % Loop through and print each channel with its max offset value
    for idx = 1:length(bad_electrodes_all_label)
        fprintf('  Channel: %-6s | Max Offset: %.2f mV\n', bad_electrodes_all_label{idx}, strange_offsets(idx));
    end
    fprintf('--------------------------------------------------\n');
else
    fprintf('All channels are within the %.1f mV threshold across all blocks.\n', threshold_mV);
end

% Log
for i_block = 1:num_block
    EEG(i_block).ALSUTRECHT.badchaninfo.offsets.offsets_mV       = offsets_mV;
    EEG(i_block).ALSUTRECHT.badchaninfo.offsets.electrodes_all   = bad_electrodes;
    EEG(i_block).ALSUTRECHT.badchaninfo.offsets.electrodes       = bad_electrodes_all_label;
    EEG(i_block).ALSUTRECHT.badchaninfo.offsets.electrodes_indx  = bad_electrodes_all_indx;
end

% -------------------------------------------------------------------------
% Adjust figure window to support a 2x2 grid
fh = figure('Name', 'Electrode Offsets Summary', 'Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', cfg.figure.visible);

% Initialize a 2x2 tiled layout
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% =========================================================================
% --- DATA PREPARATION & BAD ELECTRODE DETECTION ---
% =========================================================================
mask_plot = false(128, 1);
if ~isempty(bad_electrodes_all_indx)
    bad_electrodes_tmp = bad_electrodes_all_indx;
    bad_electrodes_tmp = bad_electrodes_tmp(bad_electrodes_tmp <= 128);
    mask_plot(bad_electrodes_tmp) = true;
end

% =========================================================================
% --- PANEL 1 (1,1): TOPOPLOT OF DC OFFSETS ---
% =========================================================================
ax1 = nexttile(tl);

% Clean up NaNs just in case to prevent topoplot crashes
data_plot = mean(offsets_mV(1:128, :), 2);
data_plot(isnan(data_plot)) = 0;
data_plot = abs(data_plot);

mytopoplot(data_plot, mask_plot, 'abs(Offsets) (mV)', ax1, threshold_mV * [0 1]);
colorbar(ax1);

% =========================================================================
% --- PANEL 2 (1,2): LONGITUDINAL BLOCK DRIFT ---
% =========================================================================
ax2 = nexttile(tl);
hold(ax2, 'on');

% Plot good channels in a subtle grey, bad channels in a distinct orange/red
for ch = 1:128
    if mask_plot(ch)
        plot(ax2, offsets_mV(ch, :), 'Color', [0.85 0.33 0.10 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        plot(ax2, offsets_mV(ch, :), 'Color', [0.2 0.2 0.2 0.25], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Visual reference boundaries for the BioSemi threshold
yline(ax2, threshold_mV, 'r--', 'LineWidth', 1.5, 'Label', 'Upper Limit');
yline(ax2, -threshold_mV, 'r--', 'LineWidth', 1.5, 'Label', 'Lower Limit');
yline(ax2, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Aesthetics for the drift panel
grid(ax2, 'on');
ax2.GridAlpha = 0.1;
ax2.Box = 'off';
xlim(ax2, [1 size(offsets_mV, 2)]);
xticks(ax2, 1:size(offsets_mV, 2));
xlabel(ax2, 'Recording Block Index');
ylabel(ax2, 'DC Offset (mV)');
title(ax2, 'Offset Trajectory Across Blocks');

if any(mask_plot)
    max_plot = max(abs(offsets_mV(1:128, :)), [], "all");
    max_plot = min(max_plot, 60);
else
    max_plot = threshold_mV;
end

ylim(max_plot * [-1 1] + [-5 5]);

% =========================================================================
% --- PANEL 3 (2,1): TOPOPLOT OF DRIFT OVER TIME ---
% =========================================================================
ax4 = nexttile(tl);

% Calculate total absolute drift (final block minus first block)
drift_per_channel = offsets_mV(1:128, end) - offsets_mV(1:128, 1);
drift_per_channel_abs = abs(drift_per_channel);
[sorted_drift, sorted_idx] = sort(drift_per_channel_abs, 'descend');

% Find a reasonable maximum for the color axis
drift_limit = max(drift_per_channel_abs);

mytopoplot(drift_per_channel, [], '\Delta Offset (mV)', ax4, drift_limit * [-1 1]);
colorbar(ax4);

% =========================================================================
% --- PANEL 4 (2,2): HIGHEST DRIFT CHANNELS ---
% =========================================================================
ax3 = nexttile(tl);

% Plot the top 10 worst drifting channels
bar(ax3, sorted_drift(1:10), 'FaceColor', [0.15 0.55 0.82], 'EdgeColor', 'none');
% set(ax3, 'XTickLabel', sorted_idx(1:10));

% Extract the channel labels for the top 10 sorted indices
chan_labels = {EEG(1).chanlocs(sorted_idx(1:10)).labels};
set(ax3, 'XTick', 1:10); % Ensure ticks match the 10 categories
set(ax3, 'XTickLabel', chan_labels);

grid(ax3, 'on');
ax3.GridAlpha = 0.1;
ax3.Box = 'off';
xlabel(ax3, 'Channel Index');
ylabel(ax3, '\Delta Offset (mV) |Block N - Block 1|');
title(ax3, 'Top 10 Channels by Total Session Drift');

% Save
save_figure(fh, EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_estimate_electrodeoffsets'], [25 20]);

fprintf('Done!\n');

end