function EEG = epoch_data(EEG,cfg)

fprintf('\n================================\n');
fprintf('Epoching data\n');
fprintf('================================\n');

if strcmpi(EEG.ALSUTRECHT.subject.task, 'MMN')
    triggers_extract = arrayfun(@(x) ['condition ' num2str(x)], cfg.mmn{1}, 'Uniformoutput', 0);
    EEG = pop_epoch(EEG, triggers_extract, cfg.mmn{2}, 'epochinfo', 'yes');
    EEG = {EEG};

elseif strcmpi(EEG.ALSUTRECHT.subject.task, 'SART')
    % EEG1: stimulus-locked (3, 6)
    % EEG2: response-locked (1)

    % Make sure that all 11s are 1s
    triggers_list = [EEG.event.edftype];
    assert(~any(triggers_list == 11));

    % SART wrt visual stimuli
    triggers_extract_1 = arrayfun(@(x) ['condition ' num2str(x)], cfg.sart1{1}, 'Uniformoutput', 0);
    EEG1 = pop_epoch(EEG, triggers_extract_1, cfg.sart1{2}, 'epochinfo', 'yes');
    EEG1.ALSUTRECHT.SART.type = 'StimulusLocked';

    % SART wrt response times
    triggers_extract_2 = arrayfun(@(x) ['condition ' num2str(x)], cfg.sart2{1}, 'Uniformoutput', 0);
    EEG2 = pop_epoch(EEG, triggers_extract_2, cfg.sart2{2}, 'epochinfo', 'yes');
    EEG2.ALSUTRECHT.SART.type  = 'ResponseLocked';

    % Organise
    EEG = {EEG1, EEG2};

elseif strcmpi(EEG.ALSUTRECHT.subject.task, 'RS') || strcmpi(EEG.ALSUTRECHT.subject.task, 'EO') || strcmpi(EEG.ALSUTRECHT.subject.task, 'EC')
    % EEG = epoch_rsdata3(EEG,cfg.rs{1},cfg.rs{2}); % OK if proc EO/EC only
    % EXT = epoch_rsdata3(EXT,cfg.rs{1},cfg.rs{2});
    % Assume one
    % EEG = epoch_rsdata2(EEG,cfg.rs{1},cfg.rs{2});   % OK if proc EO+EC together
    % EEG = {EEG};

    % Two epoching strategies:
    % 1. short (2s) epochs
    EEG1 = epoch_rsdata2(EEG, cfg.rs1{1}, cfg.rs1{2});
    % 2. long (5s) epochs
    if isfield(cfg, "rs2")
        EEG2 = epoch_rsdata2(EEG, cfg.rs2{1}, cfg.rs2{2});
    end

    % Organise
    if exist('EEG2', 'var') && ~isempty(EEG2)
        EEG = {EEG1, EEG2};
    else
        EEG = {EEG1};
    end

elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    % Maybe before ICA, as there is often a lot of noise before/after the contrations
    triggers_extract = arrayfun(@(x) ['condition ' num2str(x)], cfg.mt{1}, 'Uniformoutput', 0);
    EEG = pop_epoch(EEG, triggers_extract, cfg.mt{2}, 'epochinfo', 'yes');
    EEG = {EEG};

end

end