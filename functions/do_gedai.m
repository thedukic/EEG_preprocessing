function EEG = do_gedai(EEG, cfg)
% Reduce artifacts
% https://doi.org/10.1101/2025.10.04.680449
% https://github.com/neurotuning/GEDAI-master
% Unsupervised EEG denoising based on leadfield filtering

fprintf('\n================================\n');
fprintf('GEDAI: Unsupervised EEG denoising based on leadfield filtering\n');
fprintf('================================\n');

% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------
% artifact_threshold_type  = 'auto';
% epoch_size_in_cycles     = 12;
% lowcut_frequency         = 0.5;
% do_parallel              = true;
% do_visualise             = false;
% smoothing_window_seconds = Inf; % default

artifact_threshold_type  = cfg.gedai.artifact_threshold_type;
epoch_size_in_cycles     = cfg.gedai.epoch_size_in_cycles;
lowcut_frequency         = cfg.gedai.lowcut_frequency;
do_parallel              = cfg.gedai.do_parallel;
do_visualise             = cfg.gedai.do_visualise;
smoothing_window_seconds = cfg.gedai.smoothing_window_seconds; % default

fprintf('Smoothing window: %.1f s\n', smoothing_window_seconds);

% -------------------------------------------------------------------------
% Pre-filtering
% -------------------------------------------------------------------------
% High-pass filter at lowcut_frequency so that EEGartifacts contains no energy
% below this threshold, leaving the 0.1 - 0.5 Hz band intact upon subtraction.
EEG.data = do_ref(EEG.data, EEG.nbchan);

cfg_tmp = cfg;
cfg_tmp.flt.erp.hp = [0.5 4];
cfg_tmp.flt.rs.hp  = [0.5 4];
cfg_tmp.flt.mt.hp  = [0.5 4];

EEG_GEDAI = do_filtering(EEG, 'highpass', 'eeg', cfg_tmp);

% -------------------------------------------------------------------------
% Leadfield Preparation & Covariance Calculation
% -------------------------------------------------------------------------
load(fullfile(EEG_GEDAI.ALSUTRECHT.subject.mycodes, 'files', 'gedai', 'leadfield_5mm.mat'), "leadfield");
L = cell2mat(leadfield.leadfield(leadfield.inside).');

num_channel = size(L, 1);
L = do_ref(L, num_channel);
covRef = L * L';

% Ensure data is referenced in the same way as the leadfield
assert(ismatrix(EEG_GEDAI.data) && size(EEG_GEDAI.data, 1) == num_channel);

% Retain uncleaned state for power spectrum evaluation
EEG_GEDAI_old = EEG_GEDAI;

% -------------------------------------------------------------------------
% Run GEDAI
% -------------------------------------------------------------------------
[EEG_GEDAI_OUT, EEGartifacts, SENSAI_score, SENSAI_score_per_band, artifact_threshold_per_band, ...
    mean_ENOVA, ENOVA_per_epoch, com, ENOVA_per_band, ENOVA_per_channel] = ...
    GEDAI(EEG_GEDAI, artifact_threshold_type, epoch_size_in_cycles, lowcut_frequency, covRef, ...
    do_parallel, do_visualise, Inf, Inf, 'eeg', smoothing_window_seconds);

if all(all(round(ENOVA_per_band, 2) < 0.05))
    fprintf('Warning: Strange, GEDAI did not detect any noise. Checking the noise levels again.\n');
    artifact_threshold_type  = 'auto+';

    [EEG_GEDAI_OUT, EEGartifacts, SENSAI_score, SENSAI_score_per_band, artifact_threshold_per_band, ...
        mean_ENOVA, ENOVA_per_epoch, com, ENOVA_per_band, ENOVA_per_channel] = ...
        GEDAI(EEG_GEDAI, artifact_threshold_type, epoch_size_in_cycles, lowcut_frequency, covRef, ...
        do_parallel, do_visualise, Inf, Inf, 'eeg', smoothing_window_seconds);
end

EEG_GEDAI = EEG_GEDAI_OUT;
clearvars EEG_GEDAI_OUT

% Linearly subtract the isolated artifact signal from the main dataset.
% The 0.1 - 0.5 Hz band remains intact because EEGartifacts has no energy there.
EEG.data = EEG.data - EEGartifacts.data;

% % Check
% vis_artifacts(EEG_GEDAI, EEG_GEDAI_old);
% vis_artifacts(EEG, EEG_old);

% -------------------------------------------------------------------------
% Visualisation & Quality Control
% -------------------------------------------------------------------------
check_power_cleaning(EEG_GEDAI_old, EEG_GEDAI, [], 'gedai', cfg);

fh = plot_gedai(ENOVA_per_epoch, ENOVA_per_channel, ENOVA_per_band, artifact_threshold_type, cfg);
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_gedai_2'], [25 15]);

% -------------------------------------------------------------------------
% Logging
% -------------------------------------------------------------------------
EEG.ALSUTRECHT.gedai.SENSAI_score                 = SENSAI_score;
EEG.ALSUTRECHT.gedai.SENSAI_score_per_band        = SENSAI_score_per_band;
EEG.ALSUTRECHT.gedai.artifact_threshold_per_band  = artifact_threshold_per_band;
EEG.ALSUTRECHT.gedai.ENOVA_per_epoch              = ENOVA_per_epoch;
EEG.ALSUTRECHT.gedai.ENOVA_per_channel            = ENOVA_per_channel;
EEG.ALSUTRECHT.gedai.ENOVA_per_band               = ENOVA_per_band;
EEG.ALSUTRECHT.gedai.mean_ENOVA                   = mean_ENOVA;
EEG.ALSUTRECHT.gedai.com                          = com;

end


function data = do_ref(data, num_channel)
data_mean = sum(data, 1) / num_channel;
data = data - data_mean;
end

function fh = plot_gedai(ENOVA_per_epoch, ENOVA_per_channel, ENOVA_per_band, artifact_threshold_type, cfg)
% Initialise figure
fh = figure('Color', 'w', 'Visible', cfg.figure.visible);

% 1. ENOVA per Epoch (Temporal representation)
subplot(3, 2, [1 2]);
plot(ENOVA_per_epoch, 'LineWidth', 1.5, 'Color', [0.2, 0.4, 0.6]);
title('ENOVA per Epoch', 'FontWeight', 'bold');
xlabel('Epoch Number');
ylabel('ENOVA');
xlim([1, max(2, length(ENOVA_per_epoch))]);
grid on;
box off;

% 2. ENOVA per Channel (Spatial representation)
subplot(3, 2, [3 4]);
bar(ENOVA_per_channel, 'FaceColor', [0.4, 0.5, 0.3], 'EdgeColor', 'none');
title('ENOVA per Channel', 'FontWeight', 'bold');
xlabel('Channel Index');
ylabel('ENOVA');
xlim([0, length(ENOVA_per_channel) + 1]);
grid on;
box off;

% Topoplot (Bottom-left slot)
ax_topo = subplot(3, 2, 5);
clim_max = max(round(max(ENOVA_per_channel), 1), 0.05);
mytopoplot(ENOVA_per_channel, [], 'ENOVA per Channel', ax_topo, [0 clim_max]);
colorbar;

% Threshold label placed in the white space beneath the topoplot
% Uses normalized axis coordinates relative to the topoplot axis
text(ax_topo, 0.5, -0.22, sprintf('Threshold: %s', artifact_threshold_type), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'Color', [0.2 0.2 0.2]);

% 3. ENOVA per Band (Spectral representation)
subplot(3, 2, 6);
bar(ENOVA_per_band, 'FaceColor', [0.4, 0.2, 0.4], 'EdgeColor', 'none');
title('ENOVA per Band', 'FontWeight', 'bold');
xlabel('Frequency Band Index');
ylabel('ENOVA');
xlim([0, length(ENOVA_per_band) + 1]);
ylim_max = max(1.1 * max(ENOVA_per_band), 0.05);
ylim([0, ylim_max]);
grid on;
box off;

start_freq = 64;
num_bands  = length(ENOVA_per_band) - 1;
freqs      = start_freq * 2.^(0:-1:-(num_bands - 1));

band_names    = cell(1, length(ENOVA_per_band));
band_names{1} = 'Broadband';
for i = 1:num_bands
    band_names{i + 1} = sprintf('%g Hz', freqs(i));
end

xticks(1:length(band_names));
xticklabels(band_names);
xtickangle(45);
end