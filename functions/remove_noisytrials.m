function EEG = remove_noisytrials(EEG)

% Specify rejection parameters:
RELAX_epoching_cfg.InterpolateRejectedChannels='yes';
RELAX_epoching_cfg.SingleChannelImprobableDataThreshold=5; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
RELAX_epoching_cfg.AllChannelImprobableDataThreshold=3; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
RELAX_epoching_cfg.SingleChannelKurtosisThreshold=5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis
RELAX_epoching_cfg.AllChannelKurtosisThreshold=3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis
RELAX_epoching_cfg.reject_amp=60; % Absolute voltage amplitude threshold - if an epoch has voltages that deviate from 0 by more than this value the epoch is rejected

% Count initial epochs:
EEG.EpochRejections.InitialEpochs = size(EEG.data,3);

% Any one of these functions can be commented out to ignore those artifacts
% when creating the mask:

% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs:
ROIidx = 1:EEG.nbchan;
EEG = pop_eegthresh(EEG,1,ROIidx,-RELAX_epoching_cfg.reject_amp,RELAX_epoching_cfg.reject_amp,EEG.xmin,EEG.xmax,1,0);
EEG = pop_jointprob(EEG,1,ROIidx,RELAX_epoching_cfg.SingleChannelImprobableDataThreshold,RELAX_epoching_cfg.AllChannelImprobableDataThreshold,1,0);
EEG = pop_rejkurt(EEG,1,(1:EEG.nbchan),RELAX_epoching_cfg.SingleChannelKurtosisThreshold,RELAX_epoching_cfg.AllChannelKurtosisThreshold,1,0);
EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0);

%% If you have not filtered out data below 75Hz, you could use an objective muscle slope measure to reject epochs with remaining muscle activity:
% Use epoched data and FFT to detect slope of log frequency log
% power, add periods exceeding muscle threshold to mask. This method is
% designed for use with data that includes up to 75Hz, so  is not
% useful if frequencies below 75Hz are filtered out

% if strcmp(RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity,'yes')
%     [EEG] = RELAX_Rejecting_muscle_epochs(EEG, RELAX_epoching_cfg);
% end

end
