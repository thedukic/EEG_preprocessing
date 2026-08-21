function EEG = report_badelectrodes(EEG, cfg)

fprintf('\n================================\n');
fprintf('Generating bad electrode reports\n');
fprintf('================================\n');

% Merge all bad electrodes
EEG.ALSUTRECHT.badchaninfo.badElectrodes = unique([ ...
    EEG.ALSUTRECHT.badchaninfo.offsets.electrodes, ...
    EEG.ALSUTRECHT.badchaninfo.flat.electrodes, ...
    EEG.ALSUTRECHT.badchaninfo.prep.electrodes, ...
    EEG.ALSUTRECHT.badchaninfo.slope.electrodes]);

% -------------------------------------------------------------------------
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad  electrodes\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');

% Offsets
str = strjoin(EEG.ALSUTRECHT.badchaninfo.offsets.electrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid, 'High offset electrodes: %s\n', str);

% Flat
str = strjoin(EEG.ALSUTRECHT.badchaninfo.flat.electrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid, 'Flat electrodes: %s\n', str);

% PREP
str = strjoin(EEG.ALSUTRECHT.badchaninfo.prep.electrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid, 'PREP electrodes: %s\n', str);

% EMG
str = strjoin(EEG.ALSUTRECHT.badchaninfo.slope.electrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid, 'Shallow slope electrodes: %s\n', str);

% -------------------------------------------------------------------------
% Plot
fh = figure('Visible', cfg.figure.visible);
th = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Define
myCmap = brewermap(128, 'RdPu');
channel_mask   = strcmp({EEG.allchans.type}, 'EEG');
channel_labels = {EEG.allchans(channel_mask).labels};
channel_locs   = EEG.allchans(channel_mask);
assert(length(channel_labels) == 128);

% Create a uniform data vector matching the number of channels (all zeros)
% This ensures the map plots in a single baseline color from your colormap
num_chans = length(channel_locs);
bg_data = zeros(num_chans, 1);

% Define custom marker styling:
% Default electrodes ('emarker'): light grey small dots
% Highlighted electrodes ('emarker2'): black bold open circles (or 'x', '*', etc.)
primary_marker   = {'.', [0.7 0.7 0.7], 6, 1};
secondary_marker = {'o', 'k', 6, 2};

% --- 1. Offset ---
mask = double(ismember(channel_labels, EEG.ALSUTRECHT.badchaninfo.offsets.electrodes));
bad_indices = find(mask); % Extract the channel indices for emarker2

nexttile;
topoplot_new(bg_data, channel_locs, 'maplimits', [0 1], 'headrad', 'rim', ...
    'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'on', 'style', 'map', 'shading', 'flat', ...
    'emarker', primary_marker, 'emarker2', {bad_indices, secondary_marker{:}});
title(['Offset, N = ' num2str(sum(mask))]); axis tight;

% --- 2. Flat ---
mask = double(ismember(channel_labels, EEG.ALSUTRECHT.badchaninfo.flat.electrodes));
bad_indices = find(mask);

nexttile;
topoplot_new(bg_data, channel_locs, 'maplimits', [0 1], 'headrad', 'rim', ...
    'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'on', 'style', 'map', 'shading', 'flat', ...
    'emarker', primary_marker, 'emarker2', {bad_indices, secondary_marker{:}});
title(['Flat, N = ' num2str(sum(mask))]); axis tight;

% --- 3. PREP ---
mask = double(ismember(channel_labels, EEG.ALSUTRECHT.badchaninfo.prep.electrodes));
bad_indices = find(mask);

nexttile;
topoplot_new(bg_data, channel_locs, 'maplimits', [0 1], 'headrad', 'rim', ...
    'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'on', 'style', 'map', 'shading', 'flat', ...
    'emarker', primary_marker, 'emarker2', {bad_indices, secondary_marker{:}});
title(['PREP, N = ' num2str(sum(mask))]); axis tight;

% --- 4. EMG Slope ---
mask = double(ismember(channel_labels, EEG.ALSUTRECHT.badchaninfo.slope.electrodes));
bad_indices = find(mask);

nexttile;
topoplot_new(bg_data, channel_locs, 'maplimits', [0 1], 'headrad', 'rim', ...
    'colormap', myCmap, 'whitebk', 'on', 'electrodes', 'on', 'style', 'map', 'shading', 'flat', ...
    'emarker', primary_marker, 'emarker2', {bad_indices, secondary_marker{:}});
title(['Slope, N = ' num2str(sum(mask))]); axis tight;

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_detected_badelectrodes'], [24 7]);

fprintf('Done!\n');

end