function DATA = do_filtering(DATA, type_filter, channel_filter, cfg)

fprintf('\n================================\n');
fprintf('Filtering data (%s)\n', type_filter);
fprintf('================================\n');

% Task
task = upper(strtrim(DATA(1).ALSUTRECHT.subject.task));
fprintf('Using filter settings for the %s task.\n', task);

% Channel masks
chaneeg = find(strcmpi({DATA(1).chanlocs.type}, 'EEG'));
chanext = find(strcmpi({DATA(1).chanlocs.type}, 'EXT'));
chanemg = find(strcmpi({DATA(1).chanlocs.type}, 'EMG'));

% Logical modality flags
do_eeg = any(strcmpi(channel_filter, 'eeg'));
do_ext = any(strcmpi(channel_filter, 'ext'));
do_emg = any(strcmpi(channel_filter, 'emg'));

% -------------------------------------------------------------------------
% Clear unused edge settings to enforce filter type
% -------------------------------------------------------------------------
% Extract
cfg = cfg.flt;

switch lower(type_filter)
    case 'highpass'
        cfg.rs.lp  = [];
        cfg.mt.lp  = [];
        cfg.erp.lp = [];
        % if isfield(cfg, 'ext'), cfg.ext.lp = []; end
        cfg.emg.lp = [];

    case 'lowpass'
        cfg.rs.hp  = [];
        cfg.mt.hp  = [];
        cfg.erp.hp = [];
        % if isfield(cfg, 'ext'), cfg.ext.hp = []; end
        cfg.emg.hp = [];

    otherwise
        error('Invalid filter type: "%s". Expected ''highpass'' or ''lowpass''.', type_filter);
end

% -------------------------------------------------------------------------
% Set Task Parameters via switch-case
% -------------------------------------------------------------------------
switch task
    case {'SART', 'MMN'}
        eeg_lp = cfg.erp.lp;
        eeg_hp = cfg.erp.hp;

    case 'MT'
        eeg_lp = cfg.mt.lp;
        eeg_hp = cfg.mt.hp;

    case {'RS', 'EO', 'EC'}
        eeg_lp = cfg.rs.lp;
        eeg_hp = cfg.rs.hp;

    otherwise
        error('Unknown task: %s', task);
end

% % EXT: use dedicated cfg.ext if defined, otherwise mirror EEG task settings
% if isfield(cfg, 'ext') && (~isempty(cfg.ext.lp) || ~isempty(cfg.ext.hp))
%     ext_lp = cfg.ext.lp;
%     ext_hp = cfg.ext.hp;
% else
%     ext_lp = eeg_lp;
%     ext_hp = eeg_hp;
% end

ext_lp = eeg_lp;
ext_hp = eeg_hp;

% -------------------------------------------------------------------------
% Execute Filtering
% -------------------------------------------------------------------------
% 1. EEG
if check_filter_status('EEG', chaneeg, do_eeg, eeg_lp, eeg_hp)
    DATA = filter_signal(DATA, eeg_lp, eeg_hp, chaneeg, 'eeglab');
end

% 2. EXT
if check_filter_status('EXT', chanext, do_ext, ext_lp, ext_hp)
    DATA = filter_signal(DATA, ext_lp, ext_hp, chanext, 'eeglab');
end

% 3. EMG
if check_filter_status('EMG', chanemg, do_emg, cfg.emg.lp, cfg.emg.hp)
    DATA = filter_signal(DATA, cfg.emg.lp, cfg.emg.hp, chanemg, 'eeglab');
end

end

function should_filter = check_filter_status(type_name, chan_idx, do_flag, lp, hp)
% CHECK_FILTER_STATUS Evaluates channel availability, user request, and parameters.
% Returns true if all criteria are satisfied, otherwise prints the exact reason for skipping.

has_chans  = ~isempty(chan_idx);
has_params = ~isempty(lp) || ~isempty(hp);

if has_chans && do_flag && has_params
    should_filter = true;
    fprintf('\n--------------------------------\n');
    fprintf('%s signals:\n', upper(type_name));
    fprintf('--------------------------------\n');
else
    should_filter = false;
    if ~has_chans
        fprintf('%s filtering skipped: No matching channels found in dataset.\n', upper(type_name));
    elseif ~do_flag
        fprintf('%s filtering skipped: Turned off by user configuration.\n', upper(type_name));
    elseif ~has_params
        fprintf('%s filtering skipped: No cutoff frequencies defined.\n', upper(type_name));
    end
end

end