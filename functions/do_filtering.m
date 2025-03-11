function DATA = do_filtering(DATA,thisFiltering,cfg)

fprintf('\n================================\n');
fprintf('Filtering data (%s)\n',thisFiltering);
fprintf('================================\n');

thisTask = DATA(1).ALSUTRECHT.subject.task;
chaneeg  = find(strcmp({DATA(1).chanlocs.type},'EEG'));
chanext  = find(strcmp({DATA(1).chanlocs.type},'EXT'));
chanemg  = find(strcmp({DATA(1).chanlocs.type},'EMG'));
fprintf('Using fitler settings for the %s tasks.\n',thisTask);

% Make sure that the params are not used (just in case)
if strcmpi(thisFiltering,'highpass')
    cfg.rs.lp  = [];
    cfg.mt.lp  = [];
    cfg.erp.lp = [];
    cfg.ext.lp = [];
    cfg.emg.lp = [];

elseif strcmpi(thisFiltering,'lowpass')
    cfg.rs.hp  = [];
    cfg.mt.hp  = [];
    cfg.erp.hp = [];
    cfg.ext.hp = [];
    cfg.emg.hp = [];

else
    error('Wrong input!');
end

% EEG
if any(chaneeg)
    fprintf('EEG singals:\n');
    if strcmpi(thisTask,'SART') || strcmpi(thisTask,'MMN')
        DATA = filter_signal(DATA,cfg.erp.lp,cfg.erp.hp,chaneeg,'eeglab');
        % DATA = remove_trends(DATA);

        % EXT will be treated as EEG
        cfg.ext.lp = cfg.erp.lp;
        cfg.ext.hp = cfg.erp.hp;

    elseif strcmpi(thisTask,'MT')
        DATA = filter_signal(DATA,cfg.mt.lp,cfg.mt.hp,chaneeg,'eeglab');

        cfg.ext.lp = cfg.mt.lp;
        cfg.ext.hp = cfg.mt.hp;

    elseif strcmpi(thisTask,'RS') || strcmpi(thisTask,'EO') || strcmpi(thisTask,'EC')
        DATA = filter_signal(DATA,cfg.rs.lp,cfg.rs.hp,chaneeg,'eeglab');

        cfg.ext.lp = cfg.rs.lp;
        cfg.ext.hp = cfg.rs.hp;
    else
        error('Unknown task.');
    end
end

% EXT
if any(chanext)
    fprintf('EXT singals:\n');
    DATA = filter_signal(DATA,cfg.ext.lp,cfg.ext.hp,chanext,'eeglab');
end

% EMG
if any(chanemg)
    fprintf('EMG singals:\n');
    DATA = filter_signal(DATA,cfg.emg.lp,cfg.emg.hp,chanemg,'eeglab');
end

end