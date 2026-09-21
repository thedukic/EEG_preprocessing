function report_ica(EEG, cfg)

fprintf('\n================================\n');
fprintf('Generating ICA reports\n');
fprintf('================================\n');

% Weights for plotting
icawinv     = EEG.ALSUTRECHT.ica.icawinv;
vaf_compvar = EEG.ALSUTRECHT.ica.vaf_compvar;
num_ica     = size(icawinv, 2);

% Extract bad ICs for plotting
ICsforRemoval        = find(EEG.ALSUTRECHT.ica.final.removed);
ICsMostLikelyComplex = find(EEG.ALSUTRECHT.ica.final.complex);
ICsMostLikelyBad     = find(EEG.ALSUTRECHT.ica.final.genbad);
% ICsMostLikelyMuscle  = find(EEG.ALSUTRECHT.ica.final.muscle);
% ICsMostLikelyChannel = find(EEG.ALSUTRECHT.ica.final.channel);

% NICAtmp = length([ICsforRemoval; ICsMostLikelyMuscle; ICsMostLikelyChannel; ICsMostLikelyComplex]);
num_ica_2 = length(ICsforRemoval);

% How many rows are needed for bad ICs
% We plot the first 24 ICs
NCOL = 8;
NROW1 = 3;
num_ica_1 = NCOL * NROW1;

% Number of row for the 2nd part
NROW2 = ceil(num_ica_2 / NCOL);

% Total rows + 1 for the gap
NROW = NROW1 + NROW2 + 1;

% 1. Use square tile dimensions (e.g., 200x200 px per tile for high resolution)
tile_dim = 200;
fh = figure('Visible', cfg.figure.visible, ...
    'Position', [50, 50, NCOL * tile_dim, NROW * tile_dim], 'Color', 'w');

% 2. Use 'tight' for spacing and 'compact' for padding
tiledlayout(NROW, NCOL, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- INSIDE PART 1 LOOP ---
for i_ic = 1:min(num_ica_1, num_ica)

    % topoplot(icawinv(:, i_ic), EEG.chanlocs, ...
    %     'maplimits', max(abs(icawinv(:,i_ic)))*[-1 1], ...
    %     'headrad', 'rim', ...
    %     'colormap', myCmap1, ...
    %     'whitebk', 'on', ...
    %     'style', 'map', ...
    %     'shading', 'flat');
    mytopoplot(icawinv(:, i_ic), [], [], nexttile);

    this_label = EEG.ALSUTRECHT.ica.ICLabel.classes{EEG.ALSUTRECHT.ica.ICLabel.cvec(i_ic)};
    this_pval  = round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i_ic), 2);

    if contains(this_label, 'Brain')
        this_colour = [0.1 0.7 0.2];
    elseif contains(this_label, 'Other')
        this_colour = [0 0 0];
    else
        this_colour = [0.8 0.1 0.2];
    end

    % Reduce font size and pull title slightly closer to topoplot head rim
    t_obj = title({['#' num2str(i_ic) ' (' num2str(round(vaf_compvar(i_ic), 1)) '%)'], [this_label ', P=' num2str(this_pval)]}, ...
        'Color', this_colour, 'FontSize', 7.5, 'FontWeight', 'bold');
    t_obj.Position(2) = t_obj.Position(2) * 0.95; % Pull title slightly down toward topoplot
end

% --- INSIDE PART 2 LOOP ---
NSTART = num_ica_1 + NCOL; % Skip one full row (8 tiles) as a visual gap
cnt = 0;

for i_ic = 1:length(ICsforRemoval)
    cnt = cnt + 1;
    this_ic = ICsforRemoval(i_ic);
    this_label = EEG.ALSUTRECHT.ica.ICLabel.classes{EEG.ALSUTRECHT.ica.final.report(this_ic)};

    if any(ICsMostLikelyComplex == this_ic)
        this_label = [this_label ' (cmplx)'];
    end
    if any(ICsMostLikelyBad == this_ic)
        this_label = [this_label ' (bad)'];
    end
    if any(ICsMostLikelyBad == this_ic) && any(ICsMostLikelyComplex == this_ic)
        this_label = [this_label ' (bad&cmplx)'];
    end

    % topoplot(icawinv(:, this_ic), EEG.chanlocs, ...
    %     'maplimits', max(abs(icawinv(:, this_ic)))*[-1 1], ...
    %     'headrad', 'rim', ...
    %     'colormap', myCmap1, ...
    %     'whitebk', 'on', ...
    %     'style', 'map', ...
    %     'shading', 'flat');
    mytopoplot(icawinv(:, this_ic), [], [], nexttile(NSTART + cnt));

    t_obj = title({['#' num2str(this_ic) ' (' num2str(round(vaf_compvar(this_ic), 1)) '%)'], this_label}, ...
        'Color', [0 0 0], 'FontSize', 7.5, 'FontWeight', 'bold');
    t_obj.Position(2) = t_obj.Position(2) * 0.95;
end

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_overview'], [40 25]);

end