function EEG = do_filtering(EEG,thisFiltering,cfg)

thisTask = EEG(1).ALSUTRECHT.subject.task;
chaneeg = find(strcmp({EEG(1).chanlocs.type},'EEG'));
chanext = find(strcmp({EEG(1).chanlocs.type},'EXT'));
chanemg = find(strcmp({EEG(1).chanlocs.type},'EMG'));

if strcmpi(thisFiltering,'highpass')
    cfg.rsmt.lp = [];
    cfg.erp.lp  = [];

    % EEG
    fprintf('EEG singals:\n');
    if strcmpi(thisTask,'SART') || strcmpi(thisTask,'MMN')
        EEG = filter_signal(EEG,[],cfg.erp.hp,chaneeg,'eeglab');

    elseif strcmpi(thisTask,'MT')
        EEG = filter_signal(EEG,[],cfg.mt.hp,chaneeg,'eeglab');

    elseif strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC')
        EEG = filter_signal(EEG,[],cfg.rs.hp,chaneeg,'eeglab');

    end

    % EXT
    if any(chanext)
        fprintf('EXT singals:\n');
        EEG = filter_signal(EEG,[],cfg.ext.hp,chanext,'eeglab');
    end

    % EMG
    if strcmpi(thisTask,'MT')
        if any(chanemg)
            fprintf('EMG singals:\n');
            EEG = filter_signal(EEG,[],cfg.emg.hp,chanemg,'eeglab');
        end
    end

elseif strcmpi(thisFiltering,'lowpass')
    cfg.rsmt.hp = [];
    cfg.erp.hp  = [];

    % EEG
    fprintf('EEG singals:\n');
    if strcmpi(thisTask,'SART') || strcmpi(thisTask,'MMN')
        EEG = filter_signal(EEG,cfg.erp.lp,[],chaneeg,'eeglab');

    elseif strcmpi(thisTask,'MT')
        EEG = filter_signal(EEG,cfg.mt.lp,[],chaneeg,'eeglab');

    elseif strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC')
        EEG = filter_signal(EEG,cfg.rs.lp,[],chaneeg,'eeglab');

    end

else
    error('Wrong input!');
end

end