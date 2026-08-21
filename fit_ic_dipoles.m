function [EEG, inside_brain, good_fits] = fit_ic_dipoles(EEG, target_ics, cfg_dip)
% FIT_IC_DIPOLES
% Fits equivalent current dipoles to a specified subset of ICs and evaluates
% their location relative to the brain volume (useful for verifying EMG/muscle ICs).
%
% Inputs:
%   EEG        : EEGLAB dataset with ICA weights
%   target_ics : Array of IC indices to fit (e.g., candidate EMG components)
%   cfg_dip    : Structure containing paths, transformation settings, and plotting flags

% Thresholds
depth_threshold = -3.5; % inside brain (< 1mm)
rv_threshold = 0.15; % residual variance (< 15% == good ICs)

if nargin < 3; cfg_dip = struct(); end

% Define paths (update these if your EEGLAB installation path differs)
if ~isfield(cfg_dip, 'template_path')
    eeglab_path = fileparts(which('eeglab.m'));
    cfg_dip.template_path = fullfile(eeglab_path, 'plugins\dipfit\standard_BEM');
end
if ~isfield(cfg_dip, 'chanfile')
    cfg_dip.chanfile = 'C:\DATA\MATLAB\myCodes\preprocessing\files\elecs\biosemi128_eeglab.ced';
end

% Standard BioSemi 128 to MNI BEM transformation matrix
if ~isfield(cfg_dip, 'T')
    % cfg_dip.T = [-0.0903676 -13.1771 9.71036 0.106873 -0.00471763 -1.55223 1.0382 0.836049 0.897888];
    cfg_dip.T = [-3 -16 11, 0.00 0.00 -1.56, 1.05 0.92 0.95];
end

% Plotting flag
if ~isfield(cfg_dip, 'plot_results')
    cfg_dip.plot_results = true;
end

hdmfile = fullfile(cfg_dip.template_path, 'standard_vol.mat');
mrifile = fullfile(cfg_dip.template_path, 'standard_mri.mat');

% % Check: Launch the visualiser with the T vector applied
% chanlocs = readlocs(cfg_dip.chanfile);
% coregister(chanlocs, [], 'mesh', hdmfile, 'transform', cfg_dip.T, 'manual', 'on');

fprintf('\n=======================================================\n');
fprintf('  FITTING DIPOLES FOR %d TARGET ICs\n', length(target_ics));
fprintf('=======================================================\n');

if isempty(target_ics)
    fprintf('No target ICs provided. Skipping dipole fitting.\n');
    inside_brain = [];
    good_fits = [];
    return;
end

% 1. Initialise DIPFIT settings
EEG = pop_dipfit_settings(EEG, ...
    'hdmfile', hdmfile, ...
    'coordformat', 'MNI', ...
    'mrifile', mrifile, ...
    'chanfile', cfg_dip.chanfile, ...
    'coord_transform', cfg_dip.T, ...
    'chansel', 1:EEG.nbchan);

% 2. Fit dipoles only to the targeted ICs (to save time)
EEG = pop_multifit(EEG, target_ics, 'threshold', 100, 'dipplot', 'off', 'plotopt', {'normlen' 'on'});

% 3. Extract Residual Variance for the fitted ICs
rv_list = NaN(1, size(EEG.icaweights, 1));
for i = 1:length(target_ics)
    ic = target_ics(i);
    if isfield(EEG.dipfit.model(ic), 'rv')
        rv_list(ic) = EEG.dipfit.model(ic).rv;
    end
end

% Perform IC rejection using residual variance of the IC scalp maps.
good_fits = find(rv_list < rv_threshold);

% 4. Perform Inside/Outside Brain Criterion Check using FieldTrip
load(hdmfile, 'vol');
dipole_xyz = zeros(length(target_ics), 3);

for i = 1:length(target_ics)
    ic = target_ics(i);
    % Extract the X, Y, Z coordinates for the primary dipole of this IC
    dipole_xyz(i, :) = EEG.dipfit.model(ic).posxyz(1, :);
end

% ft_sourcedepth returns distance to the inner skull boundary
% Negative values = inside the source compartment, Positive values = outside
depth = ft_sourcedepth(dipole_xyz, vol);

% Identify components outside the brain
inside_brain = target_ics(depth <= depth_threshold);

% 5. Feedback Print
fprintf('\nDipole fitting complete.\n');
fprintf('  - Target ICs Evaluated      : %s\n', mat2str(target_ics'));
fprintf('  - Good Fits (RV < %d%%)      : %s\n', rv_threshold * 100, mat2str(good_fits));
if isempty(inside_brain)
    fprintf('  - Inside Brain ICs (Depth <= %0.1f): None\n', depth_threshold);
else
    fprintf('  - Inside Brain ICs (Depth <= %0.1f): %s\n', depth_threshold, mat2str(inside_brain));
end
fprintf('Found %d / %d targeted ICs located inside the brain volume.\n', length(inside_brain), length(target_ics));

% 6. Visualisation
if cfg_dip.plot_results
    fig_name = sprintf('Dipole Evaluation (%d ICs)', length(target_ics));
    fh = figure('Name', fig_name, 'Color', 'w');

    n_plots = length(target_ics);
    cols = ceil(sqrt(n_plots));
    rows = ceil(n_plots / cols);
    t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'IC dipole fitting & depth evaluation', 'FontSize', 14, 'FontWeight', 'bold');

    for i = 1:length(target_ics)
        ic = target_ics(i);
        ax = nexttile(t);

        % Fetch data for topoplot
        topo_data = EEG.icawinv(:, ic);
        ic_rv = rv_list(ic) * 100;
        ic_depth = depth(i);

        % Define title string and color based on depth evaluation
        if ismember(ic, inside_brain)
            loc_status = 'IN';
            title_color = [0 0.5 0]; % Green for inside
        else
            loc_status = 'OUT';
            title_color = [0.8 0 0]; % Red for outside
        end

        title_str = sprintf('IC %d | RV: %.1f%%\n%s (d=%.1f)', ic, ic_rv, loc_status, ic_depth);

        % Plot using custom topoplot function
        mytopoplot(topo_data, [], title_str, ax);

        % Apply colour to the title to easily spot rejected components
        title(ax, title_str, 'Color', title_color, 'Interpreter', 'none', 'FontSize', 10);
    end

    % 1. Determine grid layout
    n_plots = length(target_ics);
    cols    = ceil(sqrt(n_plots));
    rows    = ceil(n_plots / cols);

    % 2. Define standard dimension per topoplot tile (in cm)
    tile_w = 4.5;   % Width per tile
    tile_h = 4.5;   % Height per tile
    margin = 2.0;   % Margin for overall figure title and padding

    % 3. Calculate dynamic figure size [width, height] with sensible minimums
    fig_w = max(12, cols * tile_w);
    fig_h = max(10, rows * tile_h + margin);
    fig_dim = [fig_w, fig_h];

    % 4. Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_dipolefit'], fig_dim);
end

% Log
EEG.ALSUTRECHT.dipfit.rv_list       = rv_list;
EEG.ALSUTRECHT.dipfit.good_fits     = good_fits;
EEG.ALSUTRECHT.dipfit.depth         = depth;
EEG.ALSUTRECHT.dipfit.inside_brain = inside_brain;

end