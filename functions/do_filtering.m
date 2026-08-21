function DATA = do_filtering(DATA, type_filter, cfg)

fprintf('\n================================\n');
fprintf('Filtering data (%s)\n', type_filter);
fprintf('================================\n');

% Task
task = DATA(1).ALSUTRECHT.subject.task;
fprintf('Using filter settings for the %s task.\n', task);

% Channel masks
chaneeg  = find(strcmp({DATA(1).chanlocs.type}, 'EEG'));
chanext  = find(strcmp({DATA(1).chanlocs.type}, 'EXT'));
chanemg  = find(strcmp({DATA(1).chanlocs.type}, 'EMG'));

% Make sure that the params are not used (just in case)
if strcmpi(type_filter, 'highpass')
    cfg.rs.lp  = [];
    cfg.mt.lp  = [];
    cfg.erp.lp = [];
    cfg.ext.lp = [];
    cfg.emg.lp = [];

elseif strcmpi(type_filter, 'lowpass')
    cfg.rs.hp  = [];
    cfg.mt.hp  = [];
    cfg.erp.hp = [];
    cfg.ext.hp = [];
    cfg.emg.hp = [];

else
    error('Wrong input!');
end

% -------------------------------------------------------------------------
% EEG
if any(chaneeg)
    fprintf('\n--------------------------------\n');
    fprintf('EEG singals:\n');
    fprintf('--------------------------------\n');

    if strcmpi(task, 'SART') || strcmpi(task, 'MMN')
        DATA = filter_signal(DATA, cfg.erp.lp, cfg.erp.hp, chaneeg, 'eeglab');
        % DATA = remove_trends(DATA);

        % EXT will be treated as EEG
        cfg.ext.lp = cfg.erp.lp;
        cfg.ext.hp = cfg.erp.hp;

    elseif strcmpi(task, 'MT')
        DATA = filter_signal(DATA, cfg.mt.lp, cfg.mt.hp, chaneeg, 'eeglab');

        cfg.ext.lp = cfg.mt.lp;
        cfg.ext.hp = cfg.mt.hp;

    elseif strcmpi(task, 'RS') || strcmpi(task, 'EO') || strcmpi(task, 'EC')
        DATA = filter_signal(DATA, cfg.rs.lp, cfg.rs.hp, chaneeg, 'eeglab');

        cfg.ext.lp = cfg.rs.lp;
        cfg.ext.hp = cfg.rs.hp;
    else
        error('Unknown task.');
    end
end

% -------------------------------------------------------------------------
% EXT
if any(chanext)
    fprintf('\n--------------------------------\n');
    fprintf('EXT singals:\n');
    fprintf('--------------------------------\n');
    DATA = filter_signal(DATA, cfg.ext.lp, cfg.ext.hp, chanext, 'eeglab');
end

% -------------------------------------------------------------------------
% EMG
if any(chanemg)
    fprintf('\n--------------------------------\n');
    fprintf('EMG singals:\n');
    fprintf('--------------------------------\n');
    DATA = filter_signal(DATA, cfg.emg.lp, cfg.emg.hp, chanemg, 'eeglab');
end

end