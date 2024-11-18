function EEG = epoch_data(EEG,cfg)
% SART: EEG = [EEG1 EEG2], stimulus- and response-locked

if strcmpi(EEG.ALSUTRECHT.subject.task,'MMN')
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.mmn{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.mmn{2},'epochinfo','yes');

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    EEG0 = EEG;

    % SART wrt visual stimuli
    condLabel1 = arrayfun(@(x) ['condition ' num2str(x)],cfg.sart1{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG0,condLabel1,cfg.sart1{2},'epochinfo','yes');
    EEG.ALSUTRECHT.SART.type = 'StimulusLocked';

    % % SART wrt response times
    % condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.sart2{1},'Uniformoutput',0);
    % EEG2 = pop_epoch(EEG0,condLabel2,cfg.sart2{2},'epochinfo','yes');
    % EEG2.ALSUTRECHT.SART.type  = 'ResponseLocked';

    clearvars EEG0

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'RS') || strcmpi(EEG.ALSUTRECHT.subject.task,'EO') || strcmpi(EEG.ALSUTRECHT.subject.task,'EC')
    % EEG = epoch_rsdata3(EEG,cfg.rs{1},cfg.rs{2}); % OK if proc EO/EC only
    % EXT = epoch_rsdata3(EXT,cfg.rs{1},cfg.rs{2});
    EEG = epoch_rsdata2(EEG,cfg.rs{1},cfg.rs{2});   % OK if proc EO+EC together

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    % Maybe before ICA, as there is often a lot of noise before/after the contrations
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.mt{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.mt{2},'epochinfo','yes');

end

end