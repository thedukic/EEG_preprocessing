function [EEG, EMG, EXT] = separate_electrodetypes(DATA)

fprintf('\n================================\n');
fprintf('Separating electrode types (EEG / EMG / EXT)\n');
fprintf('================================\n');

% Ensure chanlocs and type fields exist
if ~isfield(DATA, 'chanlocs') || isempty(DATA.chanlocs) || ~isfield(DATA.chanlocs, 'type')
    error('DATA structure is missing chanlocs or chanlocs.type field.');
end

% Safely extract all channel types and labels (handling empty type entries)
all_types  = cell(1, DATA.nbchan);
all_labels = {DATA.chanlocs.labels};

for i_channel = 1:DATA.nbchan
    if isempty(DATA.chanlocs(i_channel).type)
        all_types{i_channel} = '';
    else
        all_types{i_channel} = char(DATA.chanlocs(i_channel).type);
    end
end

% 1. Separate EXT channels
ext_mask = strcmpi(all_types, 'EXT');
chanext  = all_labels(ext_mask);

if ~isempty(chanext)
    EXT  = pop_select(DATA, 'channel', chanext);
    DATA = pop_select(DATA, 'rmchannel', chanext);
    % Update local types and labels for remaining channels
    all_types(ext_mask)  = [];
    all_labels(ext_mask) = [];
else
    EXT = [];
end

% 2. Separate EMG channels
emg_mask = strcmpi(all_types, 'EMG');
chanemg  = all_labels(emg_mask);

if ~isempty(chanemg)
    EMG = pop_select(DATA, 'channel', chanemg);
    EEG = pop_select(DATA, 'rmchannel', chanemg);
else
    EMG = [];
    EEG = DATA;
end

% 3. Validate EEG set integrity
if ~isempty(EEG) && EEG.nbchan > 0
    EEG = eeg_checkset(EEG);
end
if ~isempty(EMG) && EMG.nbchan > 0
    EMG = eeg_checkset(EMG);
end
if ~isempty(EXT) && EXT.nbchan > 0
    EXT = eeg_checkset(EXT);
end

fprintf('Done! Extracted: %d EEG, %d EXT and %d EMG channels.\n', ...
    get_nbchan(EEG), get_nbchan(EXT), get_nbchan(EMG));

end

% Local helper to get channel count safely
function n = get_nbchan(EEG_struct)
if isempty(EEG_struct) || ~isfield(EEG_struct, 'nbchan')
    n = 0;
else
    n = EEG_struct.nbchan;
end
end


% function [EEG, EMG, EXT] = separate_electrodetypes(DATA)
%
% fprintf('\n================================\n');
% fprintf('Separating electrode types (EEG / EMG / EXT)\n');
% fprintf('================================\n');
%
% % Separate EXT
% chanext = {DATA.chanlocs(strcmp({DATA.chanlocs.type}, 'EXT')).labels};
% EXT = pop_select(DATA, 'channel', chanext);
% DATA = pop_select(DATA, 'rmchannel', chanext);
%
% % Separate EMG
% chanemg = {DATA.chanlocs(strcmp({DATA.chanlocs.type}, 'EMG')).labels};
% if ~isempty(chanemg)
%     EMG = pop_select(DATA, 'channel', chanemg);
%     EEG = pop_select(DATA, 'rmchannel', chanemg);
% else
%     EEG = DATA;
%     EMG = [];
% end
%
% fprintf('Done!\n');
%
% end