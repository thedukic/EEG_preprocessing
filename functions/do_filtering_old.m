function EEG = do_filtering_old(EEG,cfg)

thisTask = EEG(1).ALSUTRECHT.subject.task;
chaneeg = find(strcmp({EEG(1).chanlocs.type},'EEG'));
chanext = find(strcmp({EEG(1).chanlocs.type},'EXT'));
chanemg = find(strcmp({EEG(1).chanlocs.type},'EMG'));

% Filter EEG
if strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC') || strcmpi(thisTask,'MT')
    % RS/EO/EC/MT
    EEG = filter_signal(EEG,cfg.flt.rsmt.lp,cfg.flt.rsmt.hp,chaneeg,'eeglab');
else
    % MMN/SART (ERP)
    EEG = filter_signal(EEG,cfg.flt.erp.lp,cfg.flt.erp.hp,chaneeg,'eeglab');
end

% Filter EXT
EEG = filter_signal(EEG,cfg.flt.ext.lp,cfg.flt.ext.hp,chanext,'eeglab');

% Filter MT
if strcmpi(thisTask,'MT')
    EEG = filter_signal(EEG,cfg.flt.emg.lp,cfg.flt.emg.hp,chanemg,'eeglab');
end

end