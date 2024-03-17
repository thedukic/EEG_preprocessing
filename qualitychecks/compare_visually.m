function compare_visually(EEG,EEGRAW,theseTriggers)

% What to do with missing flat electrodes in EEGRAW?
EEG = pop_select(EEG,'nochannel',EEG.ALSUTRECHT.badchaninfo.flatElectrodes);

% Remove extreme periods
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3) && ~strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
    EEGRAW = eeg_eegrej(EEGRAW, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
end

% Epoch
if strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    EEGRAW = pop_epoch(EEGRAW,arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.sart1{1},'Uniformoutput',0),theseTriggers.sart1{2},'epochinfo','yes');
elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MMN') || strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    EEGRAW = pop_epoch(EEGRAW,arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.(lower(EEG.ALSUTRECHT.subject.task)){1},'Uniformoutput',0),theseTriggers.(lower(EEG.ALSUTRECHT.subject.task)){2},'epochinfo','yes');
elseif strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
    EEGRAW.ALSUTRECHT = EEG.ALSUTRECHT;
    EEGRAW = epoch_rsdata(EEGRAW,size(EEG.data,2)); % Epoch size 2 [s]
end
% EEGRAW = pop_rmbase(EEGRAW,[(EEGRAW.xmin)*1000 0] ,[]);

% Visual check
vis_artifacts(EEG,EEGRAW);

end