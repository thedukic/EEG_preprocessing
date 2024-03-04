function EEG = reduce_artifacts(EEG,cfgbch)
%
% Detect artifacts and try to fix them using Multi-channel Wiener filter
%
%   3 sequential MWF:
% - cleaning muscle activity first, then
% - eye blinks, then
% - horizontal eye movement and voltage drift
%

% Remove external electrodes
% EEG0 = pop_select(EEG,'nochannel',{EEG.chanlocs(~strcmp({EEG.chanlocs.type},'EEG')).labels});
% EEG0 = EEG;

% =========================================================================
% 1. Detect and correct electrode pops
% cfgbch.popType = 'large';
% EEG = detect_channelpops(EEG,cfgbch);

% =========================================================================
% 2. Detect and remove EMG
EEG = detect_channelemg(EEG,cfgbch);

% =========================================================================
% 3. Detect and remove VEOG / eye blinks
EEG = detect_eyeblinks(EEG);

% =========================================================================
% 4. Detect and remove HEOG / slow drifts
EEG = detect_channeldrifts(EEG);

% =========================================================================
% 5. Detect and remove ECG
EEG = detect_heartbeats(EEG);

% =========================================================================
% 6. Detect channel pops
% cfgbch.popType = 'small';
% EEG = detect_channelpops(EEG,cfgbch);

% =========================================================================
% Visual check
% vis_artifacts(EEG,EEG0);

end