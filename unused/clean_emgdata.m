function ALL = clean_emgdata(ALL)
%
% ALL = EEG + EXT + EMG
% 1. 50 Hz noise
% 2. MWF using ECG
% 3. ASR

% Separate EEG channels
chaneeg = {ALL.chanlocs(strcmp({ALL.chanlocs.type},'EEG')).labels};
chanext = {ALL.chanlocs(strcmp({ALL.chanlocs.type},'EXT')).labels};

EEG = pop_select(ALL,'channel',chaneeg);
EXT = pop_select(ALL,'channel',chanext);
EMG = pop_select(ALL,'rmchannel',chaneeg); % EMG+EXT

% Remove the line noise - already done?
% EMG = reduce_linenoise(EMG);

% MWF ECG
EMG = mwf_heartbeats(EMG,'EMG');

% % ASR
% EMG1 = pop_clean_rawdata(pop_select(EMG,'rmchannel',chanext),'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
%     'BurstCriterion',30,'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
% 
% EMG0 = pop_select(EMG,'rmchannel',chanext);
% vis_artifacts(EMG1,EMG0);

% Merge in the right order
EMG = pop_select(EMG,'rmchannel',chanext);
ALL = merge_eeglabsets(EEG,EXT,EMG);

end