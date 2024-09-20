function EEG = reduce_artifacts(EEG,cfgbch)
%
% Detect artifacts and minimise them using Multi-channel Wiener filter
%
% RELAXL: Three sequential MWF:
% 1. cleaning muscle activity first, then
% 2. eye blinks, then
% 3. horizontal eye movement and voltage drift

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

% A contribution made by Neil Bailey (Monash University) suggests that MWF cleaning can be improved by using sparse delay spacing. For example,
% instead of using consecutive delays [-3 -2 -1 0 1 2 3], it is now possible to specify the "delay_spacing" parameter to create a sparser sampling such as [-9 -6 -3 0 3 6 9],
% which is obtained with "delay = 3" and "delay_spacing = 3". This is also useful to include more relevant samples when your data has a higher sample rate.
%
% From informal testing by Neil Bailey, for data sampled at 1000Hz, optimal performance at cleaning eye blink artifacts was obtained by using a delay parameter setting of 8
% (so that 8 positive & negative delays are included in the MWF ), as well as a delay spacing of 16 (so that each delay is separated from the previous delay by 16 samples or 16 milliseconds).
% This provides the MWF algorithm with delay embedded covariance matrices that characterises 272ms of the data, enabling the MWF algorithm to account for a considerable proportion of each blink period.
% For muscle artifact cleaning of data sampled at 1000Hz, a delay period of 10 and delay spacing of 2 was found to be optimal.

% Good performance of the artifact removal algorithm is indicated by both high SER and high ARR.

% ASR:
% % https://academic.oup.com/sleep/article/46/12/zsad241/7275639
% Another important future task is to evaluate frequency dependency of the spatial filter in ASR: for example,
% beta and gamma band EEG power generally increases, not decreases, as a result of spatial-filter component rejection.
% This is because an output from a spatial filter is dominated by a frequency band with dominant power,
% and EEG data generally have 1/f power distribution. Thus, a spatial filter works effectively only for low-frequency signals.
% High-frequency signals are poorly decomposed and remain correlated. Rejecting correlated components causes de-cancelation
% which introduces the counterintuitive signal power increase. I reported it for the cases of ASR and independent component analysis
% -> Better then to first do ASR and try to mitigate these increases with MWF?

% % Remove external electrodes
% % EEG0 = pop_select(EEG,'nochannel',{EEG.chanlocs(~strcmp({EEG.chanlocs.type},'EEG')).labels});
% EEG0 = EEG;

% =========================================================================
% 1. MWF round 1: EMG
EEG = mwf_channelemg(EEG,cfgbch);

% =========================================================================
% 2. MWF round 2: VEOG / eye blinks
% EEG = mwf_eyeblinks(EEG);

% =========================================================================
% 3. MWF round 3: HEOG / slow drifts - hard to properly detect
% if ~strcmpi(EEG.ALSUTRECHT.subject.task,'RS') % also SART?
%     EEG = mwf_channeldrifts(EEG);
% end

% =========================================================================
% 4. MWF round 4: ECG - an overkill?
% EEG = mwf_heartbeats(EEG,EEG);

% =========================================================================
% 5. ASR - automatic, which is nice but maybe not as good as targeted MWF?
% fprintf('\nUsing ASR to fix bad segments of EEG data...\n');
%
% extchan = {EEG.chanlocs(~strcmp({EEG.chanlocs.type},'EEG')).labels};
% EXT = pop_select(EEG,'channel',extchan);
%
% originalEEG  = EEG;
% originalData = EEG.data;
% EEG = pop_clean_rawdata(pop_select(EEG,'nochannel',extchan),'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
%     'BurstCriterion',cfgbch.asr,'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
% 
% survivedDataIdx           = find(EEG.etc.clean_sample_mask);
% correspondingOriginalData = originalData(:,survivedDataIdx);
% asrPowerReductionDb       = 10*log10(var(EEG.data,0,2)./var(correspondingOriginalData,0,2));
% EEG.etc.varianceReductionInDbByAsr = asrPowerReductionDb;
% 
% EEG = merge_eeglabsets(EEG,EXT);
% fprintf('Done using ASR.\n');

% =========================================================================
% % Visual check
% EEG0.data(end,:) = 1000*EEG0.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1;
% vis_artifacts(EEG,EEG0);

end