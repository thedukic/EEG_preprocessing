function EEG = reduce_artifacts(EEG,cfgbch)
%
% Detect artifacts and minimise them using Multi-channel Wiener filter
%
% RELAXL: Three sequential MWF:
% 1. cleaning muscle activity first, then
% 2. eye blinks, then
% 3. horizontal eye movement and voltage drift
%

% The following selections determine how many epochs of each type of
% artifact to include in the mask. It's not obvious what the best choices are.
% I have selected as defaults the parameters that work best for our data.
% If only performing one MWF run, perhaps including all artifacts in the mask
% is best (as long as that leaves enough clean data for the clean data
% mask). Somers et al. (2018) suggest that it doesn't hurt to
% include clean data in the artifact mask, and so it's helpful to have wide
% boundaries around artifacts. 

% Clean periods that last for a shorter duration than the following value to be marked as artifacts, 
% and pad short artifact periods out into artifact periods of at least the following length when 
% they are shorter than this value to reduce rank deficiency issues in MWF cleaning). 
% Note that it's better to include clean periods in the artifact mask rather than the including artifact in the clean mask.


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