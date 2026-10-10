function EEG = do_asr(EEG,cfg)

fprintf('\nUsing ASR to fix bad segments of EEG data...\n');

% Separate EXT channels
eegchan = strcmp({EEG.chanlocs.type},'EEG');
extchan = {EEG.chanlocs(~eegchan).labels};
if any(extchan)
    EXT = pop_select(EEG,'channel',extchan);
    EEG = pop_select(EEG,'nochannel',extchan);
end

originalEEG  = EEG;
originalData = EEG.data;

% Can be set to run on GPU, but this option is burried in a few subfunc
% Which means that with each EEGLAB update, we need to manually change
% multiple subfunc to make it work  again ...
EEG = pop_clean_rawdata(EEG, ...
    'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
    'BurstCriterion',cfg.asr,'WindowCriterion','off','BurstRejection','off','Distance','Euclidian', ...
    'BurstCriterionRefMaxBadChns', 0, ...
    'BurstCriterionRefTolerances', [-inf 8], ...
    'WindowCriterionTolerances',[-Inf 8], 'MaxMem', 4096);

if ~isfield(EEG.etc,'clean_sample_mask')
    EEG.etc.clean_sample_mask = true(1,EEG.pnts);
end

% Check which channels have been ASR'd
originalData = originalData(eegchan,EEG.etc.clean_sample_mask);
asrPowerReductionDb = 10*log10(var(EEG.data(:,EEG.etc.clean_sample_mask),0,2) ./ var(originalData,0,2));
EEG.etc.varianceReductionInDbByAsr = asrPowerReductionDb;

fh = figure;
topoplot(asrPowerReductionDb, EEG.chanlocs(eegchan),'headrad',0.5,'whitebk','on','shading','interp','gridscale',300);
colormap(brewermap(128,'*PuOr'));
hcb = colorbar;
hcb.Title.String = "%";
title('Interpolated ASR power');

% Save
plotX=15; plotY=15;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ASR']),'-dtiff','-r300');
close(fh);

% Put back EXT channels
if any(extchan)
    EEG = merge_eeglabsets(EEG,EXT);
end

fprintf('Done using ASR.\n');

end