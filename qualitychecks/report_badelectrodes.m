function EEG = report_badelectrodes(EEG)

% Log
EEG.ALSUTRECHT.badchaninfo.badElectrodes = unique([EEG.ALSUTRECHT.badchaninfo.flatElectrodes, EEG.ALSUTRECHT.badchaninfo.PREPElectrodes, EEG.ALSUTRECHT.badchaninfo.EMGSlope]);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad  electrodes\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
str = strjoin(EEG.ALSUTRECHT.badchaninfo.flatElectrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid,'Flat electrodes: %s\n', str);
str = strjoin(EEG.ALSUTRECHT.badchaninfo.wica.fixed,', ');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA electrodes: %s\n', str);
str = arrayfun(@(x) num2str(x,'%1.1f'),EEG.ALSUTRECHT.badchaninfo.wica.pvec(EEG.ALSUTRECHT.badchaninfo.wica.ics),'uni',0);
str = strjoin(str,', ');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA P-values:   %s\n', str);
str = strjoin(EEG.ALSUTRECHT.badchaninfo.PREPElectrodes,', ');
fprintf(EEG.ALSUTRECHT.subject.fid,'PREP electrodes: %s\n', str);

% Plot
fh = figure;
th = tiledlayout(1,4);
th.TileSpacing = 'compact'; th.Padding = 'compact';

% EEG electrode colours/labels
myCmap = brewermap(128,'RdPu');
chanlabseeg = {EEG.allchans(strcmp({EEG.allchans.type},'EEG')).labels};
chanlocseeg = EEG.allchans(strcmp({EEG.allchans.type},'EEG'));

mask = double(ismember(chanlabseeg,EEG.ALSUTRECHT.badchaninfo.flatElectrodes));
nexttile;
topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map','shading','interp');
title(['Flat, N = ' num2str(sum(mask))]); axis tight;

mask = double(ismember(chanlabseeg,EEG.ALSUTRECHT.badchaninfo.wica.fixed));
nexttile;
topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map','shading','interp');
title(['wICA, N = ' num2str(sum(mask))]); axis tight;

mask = double(ismember(chanlabseeg,EEG.ALSUTRECHT.badchaninfo.PREPElectrodes));
nexttile;
topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map','shading','interp');
title(['PREP, N = ' num2str(sum(mask))]); axis tight;

mask = double(ismember(chanlabseeg,EEG.ALSUTRECHT.badchaninfo.EMGSlope));
nexttile;
topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map','shading','interp');
title(['EMG Slope, N = ' num2str(sum(mask))]); axis tight;

% Save
plotX=15; plotY=8;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_badelectrodes']),'-dtiff','-r300');
close(fh);

end