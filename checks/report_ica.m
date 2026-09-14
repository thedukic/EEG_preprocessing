function EEG = report_ica(EEG, cfg)

fprintf('\n================================\n');
fprintf('Generating ICA reports\n');
fprintf('================================\n');

% =========================================================================
% Plot the first 20 ICs + bad ICs
myCmap1 = brewermap(128, '*RdBu');
% myCmap2 = brewermap(4, 'Set1');

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

% fh = figure('Visible', cfg.figure.visible);
% th = tiledlayout(NROW, NCOL);
% th.TileSpacing = 'tight'; th.Padding = 'tight';
%
% % 1. Plot the first N ICs
% for i_ic = 1:min(num_ica_1, NICA)
%     nexttile;
%
%     % plot
%     topoplot(icawinv(:, i_ic), EEG.chanlocs, 'maplimits',max(abs(icawinv(:,i_ic)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','flat');
%
%     this_label = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i_ic)};
%     this_pval  = round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i_ic),2);
%
%     if contains(this_label, 'Brain')
%         this_colour = [0.1 0.8 0.2];
%     elseif contains(this_label, 'Other')
%         this_colour = [0 0 0];
%     else
%         this_colour = [0.8 0.1 0.2];
%     end
%
%     % Part 1 (First N ICs)
%     title({['#' num2str(i_ic) ', Var = ' num2str(round(vaf_compvar(i_ic)*100)) '%'], ...
%         [this_label ', P = ' num2str(this_pval)]}, ...
%         'Color', this_colour, 'FontSize', 8);
%
%     % axis tight;
% end
%
% % 2. Plot bad ICs
% NSTART = NICAgood + NCOL;
% cnt = 0;
% for i_type = 1:3
%     switch i_type
%         case 1
%             these_ics  = ICsforRemoval;
%             this_label = 'Removed';
%         case 2
%             these_ics  = ICsMostLikelyMuscle;
%             this_label = 'Muscle';
%         case 3
%             these_ics = ICsMostLikelyComplex;
%             this_label = 'Complex';
%     end
%     for i_ic = 1:length(these_ics)
%         cnt = cnt + 1;
%         nexttile(NSTART + cnt);
%
%         this_ic = these_ics(i_ic);
%         switch i_type
%             case 1
%                 thisLabelTmp = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.final.report(this_ic)};
%             case {2, 3}
%                 thisLabelTmp = this_label;
%         end
%
%         topoplot(icawinv(:, this_ic), EEG.chanlocs, 'maplimits',max(abs(icawinv(:,this_ic)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','flat');
%         title(['ICA' num2str(this_ic) ', ' thisLabelTmp], 'Color', myCmap2(i_type,:));
%     end
% end

% 1. Use square tile dimensions (e.g., 200x200 px per tile for high resolution)
tile_dim = 200;
fh = figure('Visible', cfg.figure.visible, ...
    'Position', [50, 50, NCOL * tile_dim, NROW * tile_dim], 'Color', 'w');

% 2. Use 'tight' for spacing and 'compact' for padding
tiledlayout(NROW, NCOL, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- INSIDE PART 1 LOOP ---
for i_ic = 1:min(num_ica_1, num_ica)
    nexttile;
    topoplot(icawinv(:, i_ic), EEG.chanlocs, ...
        'maplimits', max(abs(icawinv(:,i_ic)))*[-1 1], ...
        'headrad', 'rim', ...
        'colormap', myCmap1, ...
        'whitebk', 'on', ...
        'style', 'map', ...
        'shading', 'flat');

    this_label = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i_ic)};
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
    nexttile(NSTART + cnt);

    this_ic = ICsforRemoval(i_ic);
    this_label = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.final.report(this_ic)};

    if any(ICsMostLikelyComplex == this_ic)
        this_label = [this_label ' (cmplx)'];
    end
    if any(ICsMostLikelyBad == this_ic)
        this_label = [this_label ' (bad)'];
    end
    if any(ICsMostLikelyBad == this_ic) && any(ICsMostLikelyComplex == this_ic)
        this_label = [this_label ' (bad&cmplx)'];
    end

    topoplot(icawinv(:, this_ic), EEG.chanlocs, ...
        'maplimits', max(abs(icawinv(:, this_ic)))*[-1 1], ...
        'headrad', 'rim', ...
        'colormap', myCmap1, ...
        'whitebk', 'on', ...
        'style', 'map', ...
        'shading', 'flat');

    t_obj = title({['#' num2str(this_ic) ' (' num2str(round(vaf_compvar(this_ic), 1)) '%)'], this_label}, ...
        'Color', [0 0 0], 'FontSize', 7.5, 'FontWeight', 'bold');
    t_obj.Position(2) = t_obj.Position(2) * 0.95;
end

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_overview'], [40 25]);

end