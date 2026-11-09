function EEG = do_gedai(EEG, cfg)
% Reduce artifacts
% https://doi.org/10.1101/2025.10.04.680449
% https://github.com/neurotuning/GEDAI-master
% Unsupervised EEG denoising based on leadfield filtering

fprintf('\n================================\n');
fprintf('GEDAI: Unsupervised EEG denoising based on leadfield filtering (Two-Track 0.5 Hz)\n');
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
fprintf('Artifact threshold: %s\n', artifact_threshold_type);

% -------------------------------------------------------------------------
% Pre-filtering
% -------------------------------------------------------------------------
% High-pass filter at lowcut_frequency so that EEGartifacts contains no energy
% below this threshold, leaving the 0.1 - 0.5 Hz band intact upon subtraction.
EEG.data = do_ref(EEG.data, EEG.nbchan);

% cfg_tmp = cfg;
% cfg_tmp.flt.erp.hp = [0.5 4];
% cfg_tmp.flt.rs.hp  = [0.5 4];
% cfg_tmp.flt.mt.hp  = [0.5 4];
% EEG_GEDAI = do_filtering(EEG, 'highpass', 'eeg', cfg_tmp);

% cut_off = 0.5;
EEG_GEDAI = do_filtering_method(EEG, cfg.gedai.cut_off, cfg);

% -------------------------------------------------------------------------
% Leadfield Preparation & Covariance Calculation
% -------------------------------------------------------------------------
% Table1 = readtable('C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-10-13.txt', 'VariableNamingRule', 'preserve');
Table1 = readtable(EEG.ALSUTRECHT.subject.table1, 'VariableNamingRule', 'preserve');
Table1.ID = categorical(Table1.ID);
age = Table1.AGE(Table1.ID ==  EEG.ALSUTRECHT.subject.id);
age = mean(age);
fprintf('Average age: %.1f\n', age);

age_lf = 20:10:80;
[a, b] = min(abs(age_lf - age));
best_lf = age_lf(b);
% best_lf = 20;

file_name = ['leadfield_vrtx_' num2str(best_lf) 'yo.mat'];
file_path = fullfile('C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\template\eeg\BEM\age', file_name);
load(file_path, "leadfield");

% file_name = 'leadfield_5mm.mat';
% file_name = 'leadfield_5mm_new.mat';
% file_name = 'leadfield_vrtx.mat';
% file_path = fullfile(EEG_GEDAI.ALSUTRECHT.subject.mycodes, 'files', 'gedai', file_name);
% load(file_path, "leadfield");

leadfield = orient_leadfield(leadfield);
L = cell2mat(leadfield.leadfield(leadfield.inside).');
% L = subsample_leadfield(leadfield, 2);

indx = find(leadfield.inside, 1, "first");
num_ori = size(leadfield.leadfield{indx}, 2);
fprintf('Leadfield: %s\n', file_name);
fprintf('Using %d/%d dipoles from the leadfield.\n', size(L, 2), num_ori*sum(leadfield.inside));

num_channel = size(L, 1);
L = do_ref(L, num_channel);
covRef = L * L';

% mytopoplot(diag(covRef), []);

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
    artifact_threshold_type = 'auto+';
    fprintf('Warning: Strange, GEDAI did not detect any noise. Checking the noise levels again using %s setting.\n', artifact_threshold_type);

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

function L = subsample_leadfield(leadfield, k)
% 1. Get the original indices of the active dipoles
idx_inside = find(leadfield.inside);

% 2. Extract their 3D positions
pos_inside = leadfield.pos(idx_inside, :);

% 3. Find the unique coordinates in each dimension
x_coords = unique(pos_inside(:, 1));
y_coords = unique(pos_inside(:, 2));
z_coords = unique(pos_inside(:, 3));

% 4. Select every 2nd coordinate to convert a 5 mm grid to a 10 mm grid
% (Change the '2' to '3' for 15 mm, etc.)
x_sub = x_coords(1:k:end);
y_sub = y_coords(1:k:end);
z_sub = z_coords(1:k:end);

% 5. Create a logical mask of the dipoles that fall on this new grid
keep_mask = ismember(pos_inside(:, 1), x_sub) & ...
    ismember(pos_inside(:, 2), y_sub) & ...
    ismember(pos_inside(:, 3), z_sub);

% 6. Get the final FieldTrip indices for the subsampled dipoles
idx_subsampled = idx_inside(keep_mask);

% 7. Extract the concatenated matrix
% cell2mat automatically preserves the [Channels x 3] structure per dipole
L = cell2mat(leadfield.leadfield(idx_subsampled).');


%
x_coords_new = leadfield.pos(idx_subsampled, 1);
y_coords_new = leadfield.pos(idx_subsampled, 2);
z_coords_new = leadfield.pos(idx_subsampled, 3);
coords_new = [x_coords_new, y_coords_new, z_coords_new];

end

function lf_out = orient_leadfield(lf_in)
% ORIENT_LEADFIELD Projects unconstrained leadfields onto surface normals.
%
% Input:
%   lf_in  - FieldTrip leadfield structure with fields:
%            .leadfield {Npos x 1 cell} of [Nchan x 3]
%            .ori       {Npos x 1 cell} of [1 x 3] or [3 x 1] (unit surface normals)
%            .inside    [Npos x 1 logical]
%
% Output:
%   lf_out - Leadfield structure with oriented leadfields [Nchan x 1]

lf_out = lf_in;

inside_idx = find(lf_in.inside);

for i = 1:numel(inside_idx)
    pos_idx = inside_idx(i);

    L = lf_in.leadfield{pos_idx}; % [Nchan x 3]
    u = lf_in.ori{pos_idx};       % [1 x 3] or [3 x 1]

    if isempty(L) || isempty(u)
        continue;
    end

    % Ensure column vector [3 x 1] and unit length
    u = u(:);
    norm_u = norm(u);
    if norm_u > 0
        u = u / norm_u;
    end

    % Project: [Nchan x 3] * [3 x 1] -> [Nchan x 1]
    lf_out.leadfield{pos_idx} = L * u;
end

% Update dimord metadata
if isfield(lf_out, 'leadfielddimord')
    lf_out.leadfielddimord = '{pos}_chan';
end
end