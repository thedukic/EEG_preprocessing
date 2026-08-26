function DATA = detect_swappedelectrodes_old(DATA, cfg)

% Define
min_corr = 0.25;

% 1. Extract 3D coordinates
X = [DATA(1).chanlocs.X]';
Y = [DATA(1).chanlocs.Y]';
Z = [DATA(1).chanlocs.Z]';
channel_coords = [X, Y, Z];

% 2. Calculate Distance and Correlation Matrices
data_distance = pdist2(channel_coords, channel_coords);

% Triangulation
neighbours = find_neighbours(DATA(1).chanlocs);

% Temporary Band-pass Filter for Correlation Check
channel_eeg  = strcmp({DATA(1).chanlocs.type}, 'EEG');

num_channel = size(neighbours, 1);
assert(sum(channel_eeg) == num_channel);

EEG_TMP = filter_signal(DATA, [25 4], [5 4], 1:num_channel, 'eeglab');

% Reference
EEG_TMP = do_reref(EEG_TMP, 'aRobust');

% % Extract the data temporarily
% % data_eeg = cat(2, EEG_TMP(:).data);
% % data_eeg = data_eeg(channel_eeg, :);
%
% % Robust Spearman Rank Correlation
% % Immune to any remaining amplitude outliers
% num_block = length(EEG_TMP);
% corr_eeg = NaN(num_channel, num_channel, num_block);
% for i_block = 1:num_block
%     corr_eeg(:, :, i_block) = corr(EEG_TMP(i_block).data(channel_eeg, :)', 'Type', 'Spearman');
% end

num_block = length(EEG_TMP);
sum_z = zeros(num_channel, num_channel);
total_weight = 0;

for i_block = 1:num_block
    % 1. Get the number of samples in this specific block
    N_samples = size(EEG_TMP(i_block).data, 2);

    % 2. Calculate the block correlation
    r_block = corr(EEG_TMP(i_block).data(channel_eeg, :)', 'Type', 'Spearman');

    % 3. Convert to Fisher Z
    z_block = atanh(r_block);
    z_block(isinf(z_block)) = NaN; % Clean up the diagonal

    % 4. Define statistical weight (Degrees of freedom = N - 3)
    weight = N_samples - 3;

    % 5. Accumulate weighted Z-scores
    sum_z = sum_z + (z_block * weight);
    total_weight = total_weight + weight;
end

% 6. Divide by total weight to get the true weighted mean
mean_z = sum_z / total_weight;

% 7. Invert back to R-values
corr_eeg = tanh(mean_z);
corr_eeg(logical(eye(num_channel))) = 1; % Restore diagonal

% =========================================================================
suspected_swaps = [];
corr_eeg_local = zeros(num_channel, 1);

% 3. Automated Detection
for i_channel = 1:num_channel
    % % Select the closest 4 channels (excluding itself)
    % [sortedDistances, idx] = sort(distMatrix(i_channel,:));
    % neighbours = idx(2:5);
    % neighbourDistances = sortedDistances(2:5);
    % validNeighbours = neighbours(neighbourDistances < 40);

    validNeighbours = neighbours(i_channel, :);
    validNeighbours(validNeighbours == i_channel | validNeighbours == 0) = [];

    if length(validNeighbours) >= 2
        localCorr = corr_eeg(i_channel, validNeighbours);
        corr_eeg_local(i_channel) = mean(localCorr);

        if corr_eeg_local(i_channel) < min_corr
            suspected_swaps = [suspected_swaps, i_channel];
        end
    else
        warning('Channel %s has fewer than 2 neighbours.', DATA(1).chanlocs(i_channel).labels);
        corr_eeg_local(i_channel) = NaN;
    end
end

if ~isempty(suspected_swaps)
    fprintf('Potential swapped electrodes detected at indices: %s\n', num2str(suspected_swaps));

    % Print labels for easier identification
    labels = {DATA(1).chanlocs(suspected_swaps).labels};
    fprintf('Suspect Labels: %s\n', strjoin(labels, ', '));
else
    fprintf('No swapped electrodes detected.\n');
end

% =========================================================================
% 4. Visualisation
% =========================================================================

% Extract upper triangle to avoid duplicate pairs in the scatter plot
mask = triu(true(num_channel), 1);
dist_flat = data_distance(mask);
corr_flat = corr_eeg(mask);

% Widen and deepen the figure for a 2x4 layout
fh = figure('Name', 'Electrode Swap QA', 'Color', 'w', 'Position', [100, 100, 1400, 750], 'Visible', cfg.figure.visible);

% Setup a compact 2 row by 4 column layout
t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------------------
% ROW 1: Main Diagnostic Plots
% -------------------------------------------------------------------------

% --- Plot 1: Distance vs Correlation Scatter (Spans 1st and 2nd column) ---
ax1 = nexttile(t, [1, 2]);
hold on; box off;

scatter(dist_flat, corr_flat, 10, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);

if ~isempty(suspected_swaps)
    for s = 1:length(suspected_swaps)
        ch_idx = suspected_swaps(s);
        ch_distances = data_distance(ch_idx, setdiff(1:num_channel, ch_idx));
        ch_correlations = corr_eeg(ch_idx, setdiff(1:num_channel, ch_idx));
        scatter(ch_distances, ch_correlations, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
    end
    title('Decay Curve (Suspects Highlighted)', 'FontSize', 12, 'FontWeight', 'bold');
else
    title('Normal Correlation Decay', 'FontSize', 12, 'FontWeight', 'bold');
end

xlabel('Physical Distance (mm)');
ylabel('Signal Correlation (r)');
ylim([-1 1]);
yline(0, 'k--', 'HandleVisibility', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 10);

% --- Plot 2: Global Heatmap (Spans 3rd and 4th column) ---
ax2 = nexttile(t, [1, 2]);

topoData = corr_eeg_local;
topoData(isnan(topoData)) = 0;

mytopoplot(topoData, suspected_swaps, 'Mean Local Correlation (r)', ax2, [0.1 0.9]);
colormap(ax2, brewermap(128, 'BuPu'));
title('Local Correlation Heatmap', 'FontSize', 12, 'FontWeight', 'bold');
cb = colorbar('Location', 'EastOutside');
ylabel(cb, 'Mean Local Correlation (r)', 'FontSize', 10);

% -------------------------------------------------------------------------
% ROW 2: 4 Channels with Lowest Local Correlation
% -------------------------------------------------------------------------

% Find the 4 worst channels (ignoring actual NaNs during sorting)
[~, sorted_idx] = sort(corr_eeg_local, 'ascend', 'MissingPlacement', 'last');
worst_channels = sorted_idx(1:4);

for i = 1:4
    target_ch = worst_channels(i);
    ax_sub = nexttile(t, [1, 1]); % Each takes exactly 1 tile space

    % Extract the raw correlation vector for this specific channel
    ch_topo_data = corr_eeg(target_ch, :);

    tmp = ch_topo_data;
    tmp(target_ch) = [];
    ch_topo_data(target_ch) = max(tmp);

    % Plot individual channel correlation profile map
    mytopoplot(ch_topo_data, target_ch, '', ax_sub, [0 max(ch_topo_data)]);
    colormap(ax_sub, brewermap(128, 'YlGn')); colorbar;

    % Add label with channel name and its local correlation value
    ch_name = DATA(1).chanlocs(target_ch).labels;
    title(sprintf('%s (r = %1.2f)', ch_name, corr_eeg_local(target_ch)), ...
        'FontSize', 10, 'FontWeight', 'normal');
end

% Save
plotX = 25;
plotY = 15;
set(fh, 'InvertHardCopy', 'Off', 'Color', [1 1 1]);
set(fh, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
print(fh, fullfile(DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_detect_correlations']), '-dtiff', '-r200');
close(fh);

% -------------------------------------------------------------------------
% Log
for i_channel = 1:length(DATA)
    DATA(i_channel).ALSUTRECHT.channelcorr.corr_raw       = corr_eeg;
    DATA(i_channel).ALSUTRECHT.channelcorr.localCorrMeans = corr_eeg_local;
    DATA(i_channel).ALSUTRECHT.channelcorr.suspectedSwaps = suspected_swaps;
end

end