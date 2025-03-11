function EEG = epoch_data(EEG,cfg)

fprintf('\n================================\n');
fprintf('Epoching data\n');
fprintf('================================\n');

if strcmpi(EEG.ALSUTRECHT.subject.task,'MMN')
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.mmn{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.mmn{2},'epochinfo','yes');
    EEG = {EEG};

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    % EEG1: stimulus-locked (3, 6)
    % EEG2: response-locked (1)

    % Make sure that all 11s are 1s
    allEvents = [EEG.event.edftype];
    assert(~any(allEvents == 11));

    % Keep
    EEG0 = EEG;

    % SART wrt visual stimuli
    condLabel1 = arrayfun(@(x) ['condition ' num2str(x)],cfg.sart1{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG0,condLabel1,cfg.sart1{2},'epochinfo','yes');
    EEG.ALSUTRECHT.SART.type = 'StimulusLocked';

    % SART wrt response times
    condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.sart2{1},'Uniformoutput',0);
    EEG2 = pop_epoch(EEG0,condLabel2,cfg.sart2{2},'epochinfo','yes');
    EEG2.ALSUTRECHT.SART.type  = 'ResponseLocked';

    % Organise
    EEG = {EEG, EEG2};
    clearvars EEG0

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'RS') || strcmpi(EEG.ALSUTRECHT.subject.task,'EO') || strcmpi(EEG.ALSUTRECHT.subject.task,'EC')
    % EEG = epoch_rsdata3(EEG,cfg.rs{1},cfg.rs{2}); % OK if proc EO/EC only
    % EXT = epoch_rsdata3(EXT,cfg.rs{1},cfg.rs{2});
    EEG = epoch_rsdata2(EEG,cfg.rs{1},cfg.rs{2});   % OK if proc EO+EC together
    EEG = {EEG};

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    % Maybe before ICA, as there is often a lot of noise before/after the contrations
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.mt{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.mt{2},'epochinfo','yes');
    EEG = {EEG};

end

end