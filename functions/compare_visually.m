function compare_visually(EEG,EEGRAW,thisTask,theseTriggers)

% Remove extreme periods
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    EEGRAW = eeg_eegrej(EEGRAW, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
end

% Epoch
if strcmpi(thisTask,'SART')
    EEGRAW = pop_epoch(EEGRAW,arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.sart1{1},'Uniformoutput',0),theseTriggers.sart1{2},'epochinfo','yes');
elseif strcmpi(thisTask,'MMN') || strcmpi(thisTask,'MT')
    EEGRAW = pop_epoch(EEGRAW,arrayfun(@(x) ['condition ' num2str(x)],theseTriggers.(lower(thisTask)){1},'Uniformoutput',0),theseTriggers.(lower(thisTask)){2},'epochinfo','yes');
end
% EEGRAW = pop_rmbase(EEGRAW,[(EEGRAW.xmin)*1000 0] ,[]);

% Visual check
vis_artifacts(EEG,EEGRAW);

end