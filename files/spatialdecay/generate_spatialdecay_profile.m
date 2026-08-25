function generate_spatialdecay_profile(myPaths, subjects)

path_save = fullfile(myPaths.mycodes, 'files', 'correlation');
NSUB = length(subjects);

% Load template channel locations to initialize coordinate matrices
chanlocs = readlocs('biosemi128_eeglab.ced');
chanlbls = {chanlocs.labels};
num_channels = length(chanlbls);

% 1. Extract 3D coordinates and calculate physical distance matrix
X = [chanlocs.X]';
Y = [chanlocs.Y]';
Z = [chanlocs.Z]';
channel_coords = [X, Y, Z];
data_distance = pdist2(channel_coords, channel_coords);

% Preallocate a 3D matrix to stack correlation matrices across subjects
all_subject_corrs = NaN(num_channels, num_channels, NSUB);

fprintf('Processing %d datasets to construct the blueprint. Please wait...\n', NSUB);
numChars = 0;

for i_subj = 1:NSUB
    % Progress tracker
    fprintf(repmat('\b', 1, numChars));
    progressString = sprintf('%d/%d', i_subj, NSUB);
    fprintf('%s', progressString);
    numChars = numel(progressString);

    % File path resolution
    fileName2 = fullfile(myPaths.preproc, subjects{i_subj}, myPaths.codever, 'data', ...
        [subjects{i_subj} '_T' num2str(myPaths.visit) '_' myPaths.task '_cleandata_b.mat']);

    if exist(fileName2, 'file') == 2
        vars = whos('-file', fileName2);
        if ismember('EEG', {vars.name})
            load(fileName2, 'EEG');

            % Extract functional correlation matrix (assumes average across blocks or clean session matrix)
            subj_corr = EEG(1).ALSUTRECHT.channelcorr.corr_clean(1:num_channels, 1:num_channels);

            % Extract bad/interpolated electrode information for this specific subject
            bad_chans = EEG(1).ALSUTRECHT.badchaninfo.badElectrodes;

            % Convert bad channel labels or indices into a consistent numeric index vector
            bad_idx = [];
            if ~isempty(bad_chans)
                if iscell(bad_chans)
                    bad_idx = find(ismember({EEG(1).chanlocs.labels}, bad_chans));
                else
                    bad_idx = bad_chans;
                end
            end

            % CRITICAL MASKING STEP: Set interpolated rows/columns to NaN
            % to eliminate pre-ICA interpolation bias from the normative blueprint
            if ~isempty(bad_idx)
                subj_corr(bad_idx, :) = NaN;
                subj_corr(:, bad_idx) = NaN;
            end

            % Store into the group matrix stack
            all_subject_corrs(:, :, i_subj) = subj_corr;
        end
    end
end
fprintf('\nExtracted matrices successfully. Building group profile...\n');

% 2. Generate the unbiased aggregate blueprint using the modern omitnan syntax
group_blueprint_corr = mean(all_subject_corrs, 3, 'omitnan');

% Save the calculated blueprint matrix for direct use in the swap detection function
save(fullfile(path_save, 'normative_spatial_blueprint.mat'), 'group_blueprint_corr', 'data_distance', 'chanlocs');

% =========================================================================
% VISUALISATION: SPATIAL DECAY CURVES
% =========================================================================
% Extract the upper triangle indices to avoid self-correlation (diagonal) and duplicate pairs
upper_tri_mask = triu(true(num_channels), 1);
distance_vector = data_distance(upper_tri_mask);

fh = figure('Name', 'Spatial Decay Profile Blueprint', 'Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

% Plot individual subject profiles in thin light grey lines to visualize cohort variance
for i_subj = 1:NSUB
    subj_matrix = all_subject_corrs(:, :, i_subj);
    if ~all(isnan(subj_matrix(:)))
        subj_corr_vector = subj_matrix(upper_tri_mask);

        % Fit a local smooth curve for the subject across physical distance
        [sorted_dist, sort_idx] = sort(distance_vector);
        sorted_corr = subj_corr_vector(sort_idx);

        % Remove NaNs from individual vector before smoothing
        valid_mask = ~isnan(sorted_corr);
        if sum(valid_mask) > 10
            smoothed_corr = smooth(sorted_dist(valid_mask), sorted_corr(valid_mask), 0.2, 'loess');
            plot(sorted_dist(valid_mask), smoothed_corr, 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
        end
    end
end

% Plot the grand average blueprint profile in a bold dark blue curve
blueprint_vector = group_blueprint_corr(upper_tri_mask);
[sorted_dist, sort_idx] = sort(distance_vector);
sorted_blueprint = blueprint_vector(sort_idx);

valid_blueprint_mask = ~isnan(sorted_blueprint);
smoothed_blueprint = smooth(sorted_dist(valid_blueprint_mask), sorted_blueprint(valid_blueprint_mask), 0.1, 'loess');

plot(sorted_dist(valid_blueprint_mask), smoothed_blueprint, 'Color', [0.0 0.4 0.7], 'LineWidth', 3);

% Calculate overarching correlation between distance and grand average coupling
[overall_r, overall_p] = corr(distance_vector(valid_blueprint_mask), blueprint_vector(valid_blueprint_mask), 'Type', 'Spearman');

% Aesthetics
grid on;
ax = gca;
ax.GridAlpha = 0.15;
ax.Box = 'off';
xlabel('Physical 3D Distance between Electrodes (mm)');
ylabel('Functional Signal Correlation (Spearman r)');
title(sprintf('%s: Normative Spatial Decay Profile Baseline (N=%d)', myPaths.group, NSUB), 'FontSize', 12, 'FontWeight', 'bold');
subtitle(sprintf('Global Volume Conduction Coupling: Spearman \\rho = %.3f (p < 0.001)', overall_r), 'FontSize', 10);

% Export configuration
plotX = 20;
plotY = 15;
set(fh, 'InvertHardCopy', 'Off', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);

print(fh, fullfile(path_save, [myPaths.group '_spatial_decay_blueprint']), '-dtiff', '-r200');
close(fh);

fprintf('Final report complete. Group blueprint map saved to reports folder.\n');
end