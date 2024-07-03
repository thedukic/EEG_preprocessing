function EEG = do_filtering(EEG,thisFiltering,cfg)

thisTask = EEG(1).ALSUTRECHT.subject.task;
chaneeg = find(strcmp({EEG(1).chanlocs.type},'EEG'));
chanext = find(strcmp({EEG(1).chanlocs.type},'EXT'));
chanemg = find(strcmp({EEG(1).chanlocs.type},'EMG'));

if thisFiltering == 1
    cfg.rsmt.lp = [];
    cfg.erp.lp  = [];

    % EEG
    fprintf('EEG singals:\n');
    if strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC') || strcmpi(thisTask,'MT')
        % RS/EO/EC/MT
        EEG = filter_signal(EEG,cfg.rsmt.lp,cfg.rsmt.hp,chaneeg,'eeglab');
    else
        % MMN/SART (ERP)
        EEG = filter_signal(EEG,cfg.erp.lp,cfg.erp.hp,chaneeg,'eeglab');
    end

    % EXT
    fprintf('EXT singals:\n');
    EEG = filter_signal(EEG,cfg.ext.lp,cfg.ext.hp,chanext,'eeglab');

    % MT
    if strcmpi(thisTask,'MT')
        fprintf('EMG singals:\n');
        EEG = filter_signal(EEG,cfg.emg.lp,cfg.emg.hp,chanemg,'eeglab');
    end

elseif thisFiltering == 2
    cfg.rsmt.hp = [];
    cfg.erp.hp  = [];

    % EEG
    fprintf('EEG singals:\n');
    if strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC') || strcmpi(thisTask,'MT')
        % RS/EO/EC/MT
        EEG = filter_signal(EEG,cfg.rsmt.lp,cfg.rsmt.hp,chaneeg,'eeglab');
    else
        % MMN/SART (ERP)
        EEG = filter_signal(EEG,cfg.erp.lp,cfg.erp.hp,chaneeg,'eeglab');
    end
end

end