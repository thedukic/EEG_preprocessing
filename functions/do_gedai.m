function EEG = do_gedai(EEG, cfg)
% Reduce artifacts
% https://doi.org/10.1101/2025.10.04.680449
% https://github.com/neurotuning/GEDAI-master
% Looks at which signals are theoretically possible based on the leadfield

fprintf('\n================================\n');
fprintf('GEDAI: Unsupervised EEG denoising based on leadfield filtering\n');
fprintf('================================\n');

% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------
artifact_threshold_type = 'auto';
epoch_size_in_cycles    = 12;
lowcut_frequency        = 0.5;
do_parallel             = true;
vis_artifacts           = false;

% -------------------------------------------------------------------------
% Smoothing window
% -------------------------------------------------------------------------
smoothing_window_seconds = Inf; % default

% % Average block length
% if strcmpi(EEG.ALSUTRECHT.subject.task, 'RS')
%     smoothing_window_seconds = 3 * round(mean(EEG.ALSUTRECHT.blockinfo.block_duration));
% else
%     smoothing_window_seconds = Inf; % default
% end
fprintf('Smoothing window: %.1f s\n', smoothing_window_seconds);

% -------------------------------------------------------------------------
% Leadfield
% -------------------------------------------------------------------------
% % Leadfield with cortex does not work
% % It may be becase the cov matricesa are estiamted from a larger number of dipoles?
% load(fullfile(EEG.ALSUTRECHT.subject.mycodes, 'files', 'gedai', 'covRef_template+individual.mat'), "covRef");
% load(fullfile(EEG.ALSUTRECHT.subject.mycodes, 'files', 'gedai', 'covRef_template.mat'), "covRef");
% covRef = 'interpolated';

% Leadfield with 5 mm grid works better
load(fullfile(EEG.ALSUTRECHT.subject.mycodes, 'files', 'gedai', 'leadfield_5mm.mat'), "leadfield");
L = cell2mat(leadfield.leadfield(leadfield.inside).');

% % 1. Get the original indices of the active dipoles
% idx_inside = find(leadfield.inside);
%
% % 2. Extract their 3D positions
% pos_inside = leadfield.pos(idx_inside, :);
%
% % 3. Find the unique coordinates in each dimension
% x_coords = unique(pos_inside(:, 1));
% y_coords = unique(pos_inside(:, 2));
% z_coords = unique(pos_inside(:, 3));
%
% % 4. Select every 2nd coordinate to convert a 5 mm grid to a 10 mm grid
% % (Change the '2' to '3' for 15 mm, etc.)
% x_sub = x_coords(1:2:end);
% y_sub = y_coords(1:2:end);
% z_sub = z_coords(1:2:end);
%
% % 5. Create a logical mask of the dipoles that fall on this new grid
% keep_mask = ismember(pos_inside(:, 1), x_sub) & ...
%     ismember(pos_inside(:, 2), y_sub) & ...
%     ismember(pos_inside(:, 3), z_sub);
%
% % 6. Get the final FieldTrip indices for the subsampled dipoles
% idx_subsampled = idx_inside(keep_mask);
%
% % 7. Extract the concatenated matrix
% % cell2mat automatically preserves the [Channels x 3] structure per dipole
% L = cell2mat(leadfield.leadfield(idx_subsampled).');

num_channel = size(L, 1);
L = do_ref(L, num_channel);
covRef = L * L';

% Ensure data is ref in the same way as the leadfield
assert(ismatrix(EEG.data) && size(EEG.data, 1) == num_channel);
EEG.data = do_ref(EEG.data, num_channel);

% -------------------------------------------------------------------------
% GEDAI
% -------------------------------------------------------------------------
% GEDAI works with non-rank-deficient common-average referencing
% -> It will automatically make sure that this is the case for data (if not common-average referenced already)
% -> this is not done for leadfield covariance matrix

% Copy
EEG_old = EEG;

% Process
[EEG, EEGartifacts, SENSAI_score, SENSAI_score_per_band, artifact_threshold_per_band, mean_ENOVA, ENOVA_per_epoch, com, ENOVA_per_band, ENOVA_per_channel] = ...
    GEDAI(EEG, artifact_threshold_type, epoch_size_in_cycles, lowcut_frequency, covRef, do_parallel, vis_artifacts, ...
    Inf, Inf, 'eeg', smoothing_window_seconds);

% % Check
% vis_artifacts(EEG, EEG_old);

% -------------------------------------------------------------------------
% Visualise 0
% -------------------------------------------------------------------------
% [~, fh] = SENSAI_visualization(EEG_old, EEG, EEGartifacts, covRef, 2, 'eeg', ...
%     [], artifact_threshold_type, smoothing_window_seconds, [], cfg.figure.visible);
% SENSAI_visualization(EEG_old, EEG, EEGartifacts, covRef, 2, 'eeg', ...
%     [], artifact_threshold_type, smoothing_window_seconds);
%
% % Save
% save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_gedai_1'], [30 15]);

% -------------------------------------------------------------------------
% Visualise 1
% -------------------------------------------------------------------------
% Check
check_power_cleaning(EEG_old, EEG, [], 'gedai', cfg);

% -------------------------------------------------------------------------
% Visualise 2
% -------------------------------------------------------------------------
% Initialise figure
fh = figure('Name', 'GEDAI ENOVA Visualisation', ...
    'Color', 'w', ...
    'Position', [100, 100, 900, 800], ...
    'Visible', cfg.figure.visible);

% 1. ENOVA per Epoch (Temporal representation)
subplot(3, 2, [1 2]);
plot(ENOVA_per_epoch, 'LineWidth', 1.5, 'Color', [0.2, 0.4, 0.6]); % Dark denim blue
title('ENOVA per Epoch', 'FontWeight', 'bold');
xlabel('Epoch Number');
ylabel('ENOVA');
xlim([1, length(ENOVA_per_epoch)]);
grid on;
box off;

% 2. ENOVA per Channel (Spatial representation)
subplot(3, 2, [3 4]);
bar(ENOVA_per_channel, 'FaceColor', [0.4, 0.5, 0.3], 'EdgeColor', 'none'); % Olive green
title('ENOVA per Channel', 'FontWeight', 'bold');
xlabel('Channel Index');
ylabel('ENOVA');
xlim([0, length(ENOVA_per_channel) + 1]);
grid on;
box off;

clim_max = round(max(ENOVA_per_channel), 1);
clim_max = max(clim_max, 0.05);
mytopoplot(ENOVA_per_channel, [], 'ENOVA per Channel', subplot(3, 2, 5), [0 clim_max]);
colorbar;

% 3. ENOVA per Band (Spectral representation)
subplot(3, 2, 6);
bar(ENOVA_per_band, 'FaceColor', [0.4, 0.2, 0.4], 'EdgeColor', 'none'); % Purple
title('ENOVA per Band', 'FontWeight', 'bold');
xlabel('Frequency Band Index');
ylabel('ENOVA');
xlim([0, length(ENOVA_per_band) + 1]);
grid on;
box off;

% Dynamically generate band names based on your parameters
% nyquist = EEG.srate / 2;
% start_freq = 2^(floor(log2(nyquist)));
start_freq = 64;

% Then generate the bands:
num_bands = length(ENOVA_per_band) - 1;
freqs = start_freq * 2.^(0:-1:-(num_bands-1));

band_names = cell(1, length(ENOVA_per_band));
band_names{1} = 'Broadband';
for i = 1:num_bands
    band_names{i+1} = sprintf('%g Hz', freqs(i));
end

% Apply labels
xticks(1:length(band_names));
xticklabels(band_names);
xtickangle(45);

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_gedai_2'], [25 15]);

% -------------------------------------------------------------------------
% Log
EEG.ALSUTRECHT.gedai.SENSAI_score          = SENSAI_score;
EEG.ALSUTRECHT.gedai.SENSAI_score_per_band = SENSAI_score_per_band;
EEG.ALSUTRECHT.gedai.ENOVA_per_epoch       = ENOVA_per_epoch;
EEG.ALSUTRECHT.gedai.ENOVA_per_channel     = ENOVA_per_channel;
EEG.ALSUTRECHT.gedai.ENOVA_per_band        = ENOVA_per_band;
EEG.ALSUTRECHT.gedai.mean_ENOVA            = mean_ENOVA;
EEG.ALSUTRECHT.gedai.com                   = com;

end


function data = do_ref(data, num_channel)
% num_channel = 128;
% data_mean = sum(data, 1) / (num_channel + 1);
data_mean = sum(data, 1) / num_channel;
data = data - data_mean;
end