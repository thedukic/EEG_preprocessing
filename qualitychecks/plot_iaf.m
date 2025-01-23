function plot_iaf(paf,psdspectra1,freq1,psdspectra2,freq2,subject)
% Helper function for plotting individual alpha frequancy 

psdspectra1 = mean(psdspectra1,2);
psdspectra2 = mean(psdspectra2,2);

% Initialise
band_names   = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
band_colours = [0.9, 0.1, 0.1; 0.1, 0.9, 0.1; 0.1, 0.1, 0.9; 0.9, 0.9, 0.1; 0.9, 0.1, 0.9];
freq_bands   = make_fmask(freq1);

% Normalise for plotting
psdspectra1 = psdspectra1 ./ max(psdspectra1);
psdspectra2 = psdspectra2 ./ max(psdspectra2);

% Plot the average power spectrum
fh = figure; hold on;
plot(freq1,psdspectra1,'LineWidth',2,'Color',0.2*ones(1,3));
if ~isnan(paf)
    plot(freq2,psdspectra2,'LineWidth',1.5,'Color',0.5*ones(1,3));
end

% Plot transparent boxes for each canonical frequency band
for i = 1:length(freq_bands)
    band_mask = freq_bands{i};
    % Get the frequency indices for the band
    band_indices = find(band_mask);

    % Shade the region corresponding to the canonical band
    fill([freq1(band_indices(1)), freq1(band_indices(end)), freq1(band_indices(end)), freq1(band_indices(1))], ...
        [0, 0, max(psdspectra1), max(psdspectra1)], band_colours(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot the alpha peak frequency
if ~isnan(paf)
    [a,b] = min(abs(freq2 - paf));
    plot(freq2(b), psdspectra2(b), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
    % text(freq2(b) + 1, dataTmp2(b), sprintf('Alpha peak: %.2f Hz', paf), 'FontSize', 10, 'Color', 'k');
end

% Add labels and legend
xticks([0 4 8 13 30 50]);
xlabel('Frequency (Hz)'); ylabel('Power (a.u.)');
if ~isnan(paf)
    legend(['Average spectrum', 'IAF spectrum', band_names],'Location','eastoutside');
else
    legend(['Average spectrum', band_names],'Location','eastoutside');
end
axis tight; xlim([0 50]); pbaspect([1.618 1 1]);
title(sprintf('Alpha peak: %.2f Hz', paf));

plotX=18; plotY=8;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(subject.preproc,[subject.id '_pspectra_checks']),'-dtiff','-r300');
close(fh);

end

% =========================================================================
% Helper fuction
% =========================================================================
function freq_bands = make_fmask(freq1)
% Assume logical mask for canonical frequency bands (e.g., delta, theta, alpha, beta, gamma)
mask_delta = (freq1 > 0  & freq1 < 4);     % Delta (0-4 Hz)
mask_theta = (freq1 > 4  & freq1 < 8);     % Theta (4-8 Hz)
mask_alpha = (freq1 > 8  & freq1 < 13);    % Alpha (8-13 Hz)
mask_beta  = (freq1 > 13 & freq1 < 30);   % Beta  (13-30 Hz)
mask_gamma = (freq1 > 30 & freq1 < 48);   % Gamma (30-48 Hz)

% Combine all the masks into a cell array
freq_bands   = {mask_delta, mask_theta, mask_alpha, mask_beta, mask_gamma};
end