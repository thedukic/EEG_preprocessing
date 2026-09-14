function compare_visually(EEG, EEGRAW, theseTriggers)

% What to do with missing flat electrodes in EEGRAW?
EEG = pop_select(EEG, 'nochannel', EEG.ALSUTRECHT.badchaninfo.flatElectrodes);

% Merge blocks
if length(EEGRAW) > 1
    % EEGRAW = remove_datasetends(EEGRAW);
    EEGRAW = merge_eeglabblocks(EEGRAW);
    EEGRAW = make_extbipolar(EEGRAW);
end

% % Remove extreme periods
% if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3) && ~strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
%     EEGRAW = eeg_eegrej(EEGRAW, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
% end

% Epoch
if strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    condLabels = arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.sart1{1},'Uniformoutput',0);
    EEG    = pop_epoch(EEG,condLabels,theseTriggers.sart1{2},'epochinfo','yes');
    EEGRAW = pop_epoch(EEGRAW,condLabels,theseTriggers.sart1{2},'epochinfo','yes');

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MMN') || strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    condLabels = arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.(lower(EEG.ALSUTRECHT.subject.task)){1},'Uniformoutput',0);
    EEG    = pop_epoch(EEG,condLabels,theseTriggers.(lower(EEG.ALSUTRECHT.subject.task)){2},'epochinfo','yes');
    EEGRAW = pop_epoch(EEGRAW,condLabels,theseTriggers.(lower(EEG.ALSUTRECHT.subject.task)){2},'epochinfo','yes');

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
    % EEGRAW.ALSUTRECHT = EEG.ALSUTRECHT;
    % EEGRAW = epoch_rsdata(EEGRAW,size(EEG.data,2)); % Epoch size 2 [s]
    % EEGRAW = epoch_rsdata(EEGRAW,size(EEG.data,2)); % Epoch size 2 [s]
end
% EEGRAW = pop_rmbase(EEGRAW,[(EEGRAW.xmin)*1000 0] ,[]);

% Visual check
% EEG.etc = rmfield(EEG.etc,'clean_sample_mask');
% EEGRAW.etc = rmfield(EEGRAW.etc,'clean_sample_mask');
vis_artifacts(EEG, EEGRAW);

end