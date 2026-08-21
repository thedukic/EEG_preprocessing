function DATA = detect_swappedelectrodes(DATA, cfg)

fprintf('\n================================\n');
fprintf('Detecting swapped electrodes\n');
fprintf('================================\n');

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

EEG_TMP = filter_signal(DATA, [20 4], [3 4], 1:num_channel, 'eeglab');
EEG_TMP = do_reref(EEG_TMP, 'aRobust');

num_block = length(EEG_TMP);
sum_z = zeros(num_channel, num_channel);
total_weight = 0;

for i_block = 1:num_block
    N_samples = size(EEG_TMP(i_block).data, 2);
    r_block = corr(EEG_TMP(i_block).data(channel_eeg, :)', 'Type', 'Spearman');

    z_block = atanh(r_block);
    z_block(isinf(z_block)) = NaN;

    weight = N_samples - 3;
    sum_z = sum_z + (z_block * weight);
    total_weight = total_weight + weight;
end

mean_z = sum_z / total_weight;
corr_eeg = tanh(mean_z);
corr_eeg(logical(eye(num_channel))) = 1;

% Calculate baseline local correlation for visualisation (the blurry mess)
corr_eeg_local_baseline = zeros(num_channel, 1);
for i_channel = 1:num_channel
    validNeighbours = neighbours(i_channel, :);
    validNeighbours(validNeighbours == i_channel | validNeighbours == 0) = [];

    if length(validNeighbours) >= 2
        corr_eeg_local_baseline(i_channel) = mean(corr_eeg(i_channel, validNeighbours));
    else
        corr_eeg_local_baseline(i_channel) = NaN;
    end
end

% =========================================================================
% 3. Automated Detection: Spatial Profile Routing
% =========================================================================
fprintf('Running deterministic spatial profile routing...\n');

% Extract labels for the EEG channels to enforce intra-set boundaries
chan_labels = {DATA(1).chanlocs(1:num_channel).labels};

% OLD
% [resolved_swaps, corr_eeg_corrected] = resolve_all_swaps_profile(corr_eeg, data_distance, chan_labels);

% Load the normative blueprint matrix (adjust path as needed)
blueprint_file = fullfile(DATA(1).ALSUTRECHT.subject.mycodes, 'files', 'correlation', 'normative_spatial_blueprint.mat');
if exist(blueprint_file, 'file')
    load(blueprint_file, 'group_blueprint_corr');
else
    error('Normative blueprint not found. Run generate_spatialdecay_profile first.');
end

% Pass the blueprint INSTEAD of the data_distance matrix
[resolved_swaps, corr_eeg_corrected] = resolve_all_swaps_profile(corr_eeg, group_blueprint_corr, chan_labels);
flat_swaps = resolved_swaps(:)';

if ~isempty(resolved_swaps)
    for s = 1:size(resolved_swaps, 1)
        chA = resolved_swaps(s, 1);
        chB = resolved_swaps(s, 2);

        fprintf('CONFIRMED SWAP: Automatically correcting metadata for %s and %s\n', ...
            DATA(1).chanlocs(chA).labels, DATA(1).chanlocs(chB).labels);

        % Apply the fix to the physical locations across all blocks
        for i_blk = 1:length(DATA)
            temp_loc = DATA(i_blk).chanlocs(chA);

            % Move B into A
            DATA(i_blk).chanlocs(chA).X = DATA(i_blk).chanlocs(chB).X;
            DATA(i_blk).chanlocs(chA).Y = DATA(i_blk).chanlocs(chB).Y;
            DATA(i_blk).chanlocs(chA).Z = DATA(i_blk).chanlocs(chB).Z;
            DATA(i_blk).chanlocs(chA).labels = DATA(i_blk).chanlocs(chB).labels;
            if isfield(DATA(i_blk).chanlocs, 'theta')
                DATA(i_blk).chanlocs(chA).theta = DATA(i_blk).chanlocs(chB).theta;
                DATA(i_blk).chanlocs(chA).radius = DATA(i_blk).chanlocs(chB).radius;
                DATA(i_blk).chanlocs(chA).sph_theta = DATA(i_blk).chanlocs(chB).sph_theta;
                DATA(i_blk).chanlocs(chA).sph_phi = DATA(i_blk).chanlocs(chB).sph_phi;
                DATA(i_blk).chanlocs(chA).sph_radius = DATA(i_blk).chanlocs(chB).sph_radius;
            end

            % Move original A into B
            DATA(i_blk).chanlocs(chB).X = temp_loc.X;
            DATA(i_blk).chanlocs(chB).Y = temp_loc.Y;
            DATA(i_blk).chanlocs(chB).Z = temp_loc.Z;
            DATA(i_blk).chanlocs(chB).labels = temp_loc.labels;
            if isfield(DATA(i_blk).chanlocs, 'theta')
                DATA(i_blk).chanlocs(chB).theta = temp_loc.theta;
                DATA(i_blk).chanlocs(chB).radius = temp_loc.radius;
                DATA(i_blk).chanlocs(chB).sph_theta = temp_loc.sph_theta;
                DATA(i_blk).chanlocs(chB).sph_phi = temp_loc.sph_phi;
                DATA(i_blk).chanlocs(chB).sph_radius = temp_loc.sph_radius;
            end
        end
    end
else
    fprintf('No swapped electrodes detected. Array is structurally sound.\n');
end

% =========================================================================
% 4. Visualisation (Before & After QA Report)
% =========================================================================

% 1. Calculate the spatial match (\rho) before and after correction
rho_before = zeros(num_channel, 1);
rho_after = zeros(num_channel, 1);

for i = 1:num_channel
    dists = data_distance(i, :);

    % Before
    corrs_b = corr_eeg(i, :);
    idx_b = setdiff(1:num_channel, i);
    rho_before(i) = corr(dists(idx_b)', corrs_b(idx_b)', 'Type', 'Spearman', 'Rows', 'complete');

    % After (If no swaps, this is identical to before)
    if exist('corr_eeg_corrected', 'var')
        corrs_a = corr_eeg_corrected(i, :);
        rho_after(i) = corr(dists(idx_b)', corrs_a(idx_b)', 'Type', 'Spearman', 'Rows', 'complete');
    else
        rho_after(i) = rho_before(i);
        corr_eeg_corrected = corr_eeg;
    end
end

% Set up the 2x2 layout
fh = figure('Name', 'Electrode Swap QA', 'Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', cfg.figure.visible);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Plot 1: Array Health (Before) ---
d_min = round(max(rho_before), 1);
d_max = round(min(rho_before), 1);

ax1 = nexttile(t);
mytopoplot(rho_before, flat_swaps, 'Spatial Decay (\rho)', ax1, [d_max d_min]);
colormap(ax1, brewermap(128, '*BuPu'));
title('Array Health BEFORE Correction', 'FontSize', 12, 'FontWeight', 'bold');
cb1 = colorbar; ylabel(cb1, 'Spearman \rho');

% --- Plot 2: Array Health (After) ---
ax2 = nexttile(t);
mytopoplot(rho_after, [], 'Spatial Decay (\rho)', ax2, [d_max d_min]);
colormap(ax2, brewermap(128, '*BuPu'));
title('Array Health AFTER Correction', 'FontSize', 12, 'FontWeight', 'bold');
cb2 = colorbar; ylabel(cb2, 'Spearman \rho');

% --- Plot 3 & 4: The Smoking Gun Footprint ---
ax3 = nexttile(t);
ax4 = nexttile(t);

if ~isempty(flat_swaps)
    % Pick the first swapped channel to demonstrate the fix
    target_ch = flat_swaps(1);
    ch_name = DATA(1).chanlocs(target_ch).labels;

    % Pre-Fix Footprint
    ch_topo_before = corr_eeg(target_ch, :);
    tmp = ch_topo_before; tmp(target_ch) = []; ch_topo_before(target_ch) = max(tmp);

    mytopoplot(ch_topo_before, target_ch, '', ax3, [0 max(ch_topo_before)]);
    colormap(ax3, brewermap(128, 'YlGn'));
    title(sprintf('Signal Footprint BEFORE: %s', ch_name), 'FontSize', 12, 'FontWeight', 'bold');

    % Post-Fix Footprint
    ch_topo_after = corr_eeg_corrected(target_ch, :);
    tmp = ch_topo_after; tmp(target_ch) = []; ch_topo_after(target_ch) = max(tmp);

    mytopoplot(ch_topo_after, target_ch, '', ax4, [0 max(ch_topo_after)]);
    colormap(ax4, brewermap(128, 'YlGn'));
    title(sprintf('Signal Footprint AFTER: %s', ch_name), 'FontSize', 12, 'FontWeight', 'bold');
else
    % Fallback if the array is perfectly healthy
    % Clear and hide the axes completely without forcing visibility
    axis(ax3, 'off');
    title(ax3, 'No Swaps Detected', 'FontSize', 12, 'FontWeight', 'bold');

    axis(ax4, 'off');
    title(ax4, 'Array is Structurally Sound', 'FontSize', 12, 'FontWeight', 'bold');
end

% Save and log
% plotX = 25; plotY = 20;
% set(fh, 'InvertHardCopy', 'Off', 'Color', [1 1 1], 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY]);
% print(fh, fullfile(DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_detect_correlations']), '-dtiff', '-r200'); close(fh);
save_figure(fh, DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_detect_correlations'], [25 20]);

% -------------------------------------------------------------------------
for i_channel = 1:length(DATA)
    DATA(i_channel).ALSUTRECHT.channelcorr.corr_raw          = corr_eeg;
    DATA(i_channel).ALSUTRECHT.channelcorr.corr_corrected    = corr_eeg_corrected;
    DATA(i_channel).ALSUTRECHT.channelcorr.localCorrBaseline = corr_eeg_local_baseline;
    DATA(i_channel).ALSUTRECHT.channelcorr.resolvedSwaps     = resolved_swaps;
end

end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function [swaps, current_corr] = resolve_all_swaps_profile(corr_matrix, blueprint_corr, labels)
% Resolves swaps by evaluating how perfectly a signal's correlation profile
% aligns with the normative blueprint profile for that specific scalp location.
N = size(corr_matrix, 1);
swaps = [];
current_corr = corr_matrix;

% Extract the alphabetical prefix (e.g., 'A', 'B', 'C', 'D') from Biosemi labels
prefixes = cell(N, 1);
for k = 1:N
    prefixes{k} = regexp(labels{k}, '^[A-Za-z]+', 'match', 'once');
end

% Precompute valid intra-set pairs to optimize the loop
valid_pair = false(N, N);
for i = 1:N
    for j = (i+1):N
        valid_pair(i, j) = strcmp(prefixes{i}, prefixes{j}) && ~isempty(prefixes{i});
    end
end

while true
    Benefit = zeros(N, N);

    for i = 1:N
        for j = (i+1):N
            if ~valid_pair(i, j)
                continue;
            end

            % Exclude the two targets being tested so they don't bias the fit
            idx = setdiff(1:N, [i, j]);

            % Baseline: How well do the signals match the blueprint for their CURRENT locations?
            % (Positive Spearman: because a perfect functional match = 1.0)
            match_ii = corr(current_corr(i, idx)', blueprint_corr(i, idx)', 'Type', 'Spearman', 'Rows', 'complete');
            match_jj = corr(current_corr(j, idx)', blueprint_corr(j, idx)', 'Type', 'Spearman', 'Rows', 'complete');

            % Hypothetical: How well would they match the blueprints if we SWAPPED them?
            match_ij = corr(current_corr(i, idx)', blueprint_corr(j, idx)', 'Type', 'Spearman', 'Rows', 'complete');
            match_ji = corr(current_corr(j, idx)', blueprint_corr(i, idx)', 'Type', 'Spearman', 'Rows', 'complete');

            % Net mathematical benefit of performing this swap
            Benefit(i, j) = match_ij + match_ji - match_ii - match_jj;
        end
    end

    [max_benefit, max_idx] = max(Benefit(:));

    if max_benefit > 0.35
        [chA, chB] = ind2sub([N, N], max_idx);
        swaps = [swaps; chA, chB];
        fprintf('Profile Routing identified swap between index %d and %d (Net Rank Benefit: +%.3f)\n', chA, chB, max_benefit);

        % Virtually swap the correlation matrix rows/cols to allow detection of multiple swaps
        temp_corr = current_corr;
        temp_corr([chA, chB], :) = temp_corr([chB, chA], :);
        temp_corr(:, [chA, chB]) = temp_corr(:, [chB, chA]);
        current_corr = temp_corr;
    else
        break;
    end
end
end

% function [swaps, current_corr] = resolve_all_swaps_profile(corr_matrix, data_distance, labels)
% % Resolves swaps by evaluating how perfectly a signal's correlation profile
% % aligns with a location's distance profile across the ENTIRE array.
% % Restricts hypothetical swaps to electrodes within the same prefix set.
% N = size(corr_matrix, 1);
% swaps = [];
% current_corr = corr_matrix;
%
% % Extract the alphabetical prefix (e.g., 'A', 'B', 'C', 'D') from Biosemi labels
% prefixes = cell(N, 1);
% for k = 1:N
%     prefixes{k} = regexp(labels{k}, '^[A-Za-z]+', 'match', 'once');
% end
%
% % Precompute valid intra-set pairs to optimize the loop
% valid_pair = false(N, N);
% for i = 1:N
%     for j = (i+1):N
%         % Only allow pairs that share the exact same alphabetical prefix
%         valid_pair(i, j) = strcmp(prefixes{i}, prefixes{j}) && ~isempty(prefixes{i});
%     end
% end
%
% while true
%     Benefit = zeros(N, N);
%
%     for i = 1:N
%         for j = (i+1):N
%             % Skip pairs that belong to different Biosemi sets
%             if ~valid_pair(i, j)
%                 continue;
%             end
%
%             % Exclude the two targets being tested so they don't bias the fit
%             idx = setdiff(1:N, [i, j]);
%
%             % Baseline: How well do the signals currently match their physical coordinates?
%             % (Negative Spearman: because as distance goes UP, correlation goes DOWN)
%             match_ii = -corr(current_corr(i, idx)', data_distance(i, idx)', 'Type', 'Spearman', 'Rows', 'complete');
%             match_jj = -corr(current_corr(j, idx)', data_distance(j, idx)', 'Type', 'Spearman', 'Rows', 'complete');
%
%             % Hypothetical: How well would they match if we swapped them?
%             match_ij = -corr(current_corr(i, idx)', data_distance(j, idx)', 'Type', 'Spearman', 'Rows', 'complete');
%             match_ji = -corr(current_corr(j, idx)', data_distance(i, idx)', 'Type', 'Spearman', 'Rows', 'complete');
%
%             % Net mathematical benefit of performing this swap
%             Benefit(i, j) = match_ij + match_ji - match_ii - match_jj;
%         end
%     end
%
%     [max_benefit, max_idx] = max(Benefit(:));
%
%     % A gain of > 0.05 across 126 nodes is a mathematical certainty of a structural fix.
%     if max_benefit > 0.05
%         [chA, chB] = ind2sub([N, N], max_idx);
%         swaps = [swaps; chA, chB];
%         fprintf('Profile Routing identified swap between index %d and %d (Net Rank Benefit: +%.3f)\n', chA, chB, max_benefit);
%
%         % Virtually swap the correlation matrix rows/cols to allow detection of multiple swaps
%         temp_corr = current_corr;
%         temp_corr([chA, chB], :) = temp_corr([chB, chA], :);
%         temp_corr(:, [chA, chB]) = temp_corr(:, [chB, chA]);
%         current_corr = temp_corr;
%     else
%         break; % No further structural improvements found
%     end
% end
% end