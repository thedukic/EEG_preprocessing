function report_badsegments(EEG, mask_noise, tag_figure, opt_figure)
% REPORT_BADSEGMENTS Visualise and export excluded data intervals across blocks
%
% Inputs:
%   EEG        : EEGLAB dataset struct or array of structs EEG(1:N)
%   mask_noise : Cell array of size [1 x num_block], where each cell contains:
%                - an [N x 2] matrix of [start_sample, stop_sample], OR
%                - a [1 x n_samples] logical mask
%   tag_figure : String identifier for the saved figure (default: 'bad_segments')
%   opt_figure : Figure visibility ('on' or 'off', default: 'on')

if nargin < 3 || isempty(tag_figure), tag_figure = 'bad_segments'; end
if nargin < 4 || isempty(opt_figure), opt_figure = 'on'; end

num_block = length(EEG);
assert(length(mask_noise) == num_block, ...
    'Length of mask_noise (%d) must match number of EEG datasets (%d).', ...
    length(mask_noise), num_block);

% Calculate global session summary
total_samples_all = sum([EEG.pnts]);
rej_samples_all   = 0;

% Create figure container
fh = figure('Visible', opt_figure, 'Color', 'w');
th = tiledlayout(num_block, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

for i = 1:num_block
    srate = EEG(i).srate;
    pnts  = EEG(i).pnts;
    time_factor = srate * 60; % Convert samples to minutes
    total_min   = pnts / time_factor;

    % 1. Standardise intervals input (accepts [N x 2] or 1D logical mask)
    block_input = mask_noise{i};
    if islogical(block_input) || (isvector(block_input) && numel(block_input) == pnts)
        diff_mask = diff([0, block_input(:)', 0]);
        intervals = [find(diff_mask == 1)', find(diff_mask == -1)' - 1];
    else
        intervals = block_input; % Assumed [N x 2]
    end

    % 2. Calculate rejection metrics for block
    if ~isempty(intervals)
        n_bursts    = size(intervals, 1);
        rej_samples = sum(intervals(:, 2) - intervals(:, 1) + 1);
    else
        n_bursts    = 0;
        rej_samples = 0;
    end
    pct_rej = (rej_samples / pnts) * 100;
    rej_samples_all = rej_samples_all + rej_samples;

    % 3. Render block track
    ax = nexttile;
    hold(ax, 'on');

    % Background bar (Clean / Retained Data)
    patch(ax, [0 total_min total_min 0], [0 0 1 1], ...
        [0.70 0.82 0.95], 'EdgeColor', 'none');

    % Foreground patches (Rejected Data, Vectorised)
    if n_bursts > 0
        x_starts = intervals(:, 1)' / time_factor;
        x_stops  = intervals(:, 2)' / time_factor;

        % Build coordinate matrices: 4 vertices per rectangular segment
        X_coords = [x_starts; x_stops; x_stops; x_starts];
        Y_coords = repmat([0; 0; 1; 1], 1, n_bursts);

        patch(ax, X_coords, Y_coords, [0.90 0.25 0.25], ...
            'FaceAlpha', 0.65, 'EdgeColor', 'none');
    end

    % 4. Format Axes
    set(ax, 'YLim', [0 1], 'YTick', [], 'XLim', [0 total_min]);
    box(ax, 'on');
    grid(ax, 'on');

    % Informative Y-label per tile
    block_title = sprintf('Block %d', i);
    ylabel(ax, sprintf('%s\n[%.1f%% rej]', block_title, pct_rej), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 9);

    if i == num_block
        xlabel(ax, 'Time (minutes)', 'FontWeight', 'bold');
    else
        set(ax, 'XTickLabel', []);
    end

    hold(ax, 'off');
end

% Global Title with session totals
global_pct_rej = (rej_samples_all / total_samples_all) * 100;
title(th, sprintf('Data Exclusion Summary: %.1f%% Total Data Excluded', global_pct_rej), 'FontWeight', 'bold', 'FontSize', 12);

% Save output if metadata exists
fig_dir = EEG(1).ALSUTRECHT.subject.figures;
sub_id  = EEG(1).ALSUTRECHT.subject.id;
save_figure(fh, fig_dir, [sub_id '_detected_' tag_figure], [15, max(4, num_block * 2.2)]);

% Clean up headless handles to prevent memory leaks in batch mode
if strcmpi(opt_figure, 'off')
    close(fh);
end

end