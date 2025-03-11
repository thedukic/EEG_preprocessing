function [EEG, NumberTrials] = detect_badepochs(EEG,cfg)

% =========================================================================
% Data quaility checks
% =========================================================================

fprintf('\n================================\n');
fprintf('Detecting bad epochs\n');
fprintf('================================\n');

% Note the number of trials
NumberTrials = NaN(4,1);
NumberTrials(1) = size(EEG.data,3);

% =============================================
% A. EEGLAB-based rejection
% Any one of these functions can be commented out to ignore those artifacts when creating the mask
% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs

% ROIidx = 1:128; % Use only EEG electrodes!
% fprintf('\n--------------------------------\n');
% fprintf('Max. amplitude (>abs(%d uV))\n',cfg.epoch.rejectAmp);
% fprintf('--------------------------------\n');
% EEG = pop_eegthresh(EEG,1,ROIidx,-cfg.epoch.rejectAmp,cfg.epoch.rejectAmp,EEG.xmin,EEG.xmax,1,0);
%
% fprintf('\n--------------------------------\n');
% fprintf('Improbable data\n');
% fprintf('--------------------------------\n');
% EEG = pop_jointprob(EEG,1,ROIidx,cfg.epoch.singleChannelImprobableDataThreshold,cfg.epoch.allChannelImprobableDataThreshold,1,0);
%
% fprintf('\n--------------------------------\n');
% fprintf('Kurtosis\n');
% fprintf('--------------------------------\n');
% EEG = pop_rejkurt(EEG,1,ROIidx,cfg.epoch.singleChannelKurtosisThreshold,cfg.epoch.allChannelKurtosisThreshold,1,0);
%
% fprintf('\n--------------------------------\n');
% fprintf('Combining and rejecting\n');
% fprintf('--------------------------------\n');
% EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
% EEG = pop_rejepoch(EEG, EEG.reject.rejglobal, 0);

fprintf('\n--------------------------------\n');
fprintf('Max. amplitude (>abs(%d uV))\n',cfg.epoch.rejectAmp);
fprintf('--------------------------------\n');
% It could be that some people have very strong oscillations
% Example: Resting-state alpha oscillations
% We do not want to exclude these data
% -> bandstop filter 5-25 Hz
% -> EOG <5 Hz
% -> EMG >25 Hz
freqStop = [2 20]; % Hz

fprintf('Temporarily bandstop filtering [%d-%d Hz] the data.\n',freqStop);
fprintf('This prevents removal of trials with strong brain oscillations.\n');
eegchans = strcmp({EEG.chanlocs.type},'EEG');
dataTmp  = double(EEG.data(eegchans,:));
[b, a]   = butter(5, freqStop / (EEG.srate / 2), 'stop');
assert(isstable(b,a), 'Bandstop filter unstable.');
dataTmp  = filtfilt(b, a, dataTmp')';
dataTmp  = dataTmp - mean(dataTmp,1);
dataTmp  = dataTmp - mean(dataTmp,2);
dataTmp  = reshape(dataTmp, size(EEG.data(eegchans,:,:)));
badTrialTreshold = squeeze(any(abs(dataTmp) > cfg.epoch.rejectAmp, [1 2]));

% badTrialTmp = find(badTrialTreshold);
% figure; tiledlayout(1,2);
% mytopoplot(mean(dataTmp(:,:,badTrialTmp).^2,[2 3]),[],'Filtered',nexttile); colorbar;
% dataTmp = double(EEG.data(eegchans,:,:));
% mytopoplot(mean(dataTmp(:,:,badTrialTmp).^2,[2 3]),[],'Raw',nexttile); colorbar;
% disp(sum(badTrialTreshold));
%
% [NCHN,NPTS,NTRL]= size(dataTmp);
% for i = 1:NTRL
%     [psdspectra(:,:,i), freq] = pwelch(dataTmp(:,:,i)',NPTS,0,NPTS,EEG.srate);
% end
% figure; plot(freq, mean(psdspectra,3));

if any(badTrialTreshold)
    EEG = pop_rejepoch(EEG,badTrialTreshold,0);
else
    fprintf('No high voltage trials are found.\n');
end

% Note the number of trials
NumberTrials(2) = EEG.trials;

% =============================================
% B. EMG-slope-based rejection
fprintf('\n--------------------------------\n');
fprintf('EMG slopes\n');
fprintf('--------------------------------\n');

% % 0.
% slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
% slopesChannelsxEpochs = slopesChannelsxEpochs > cfg.bch.muscleSlopeThreshold;
% BadEpochs = sum(slopesChannelsxEpochs, 1);
% if mean(BadEpochs) > 0.05
%     EEG = denoise_emg(EEG);
% end

% 1. Interpolate
slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
slopesChannelsxEpochs = slopesChannelsxEpochs > cfg.bch.muscleSlopeThreshold;

% Interpolate trials that do not have a lot of electrodes contaminated by EMG
if any(slopesChannelsxEpochs,"all")
    % Organise input for interpolation
    badElecsPerTrial = arrayfun(@(x) find(slopesChannelsxEpochs(:,x)), 1:EEG.trials, 'UniformOutput', false);
    eegchans = strcmp({EEG.chanlocs.type},'EEG');
    chanLocs = EEG.chanlocs(eegchans);

    % Max number of contaminated electrodes
    fprintf('The maximum number of EMG-contaminated electrodes: %d\n', cfg.epoch.NinterpMax);
    [data, report] = interpolate_epochs(EEG.data(eegchans,:,:),chanLocs,badElecsPerTrial,[],cfg.epoch.NinterpMax);

    % EEGNEW = EEG;
    % EEGNEW.data(eegchans,:,:) = data;
    % vis_artifacts(EEGNEW,EEG);

    % Return the data to the struct
    EEG.data(eegchans,:,:) = data;

    % These are too noisy to be saved
    badTrialMuscle = false(1,EEG.trials);
    badTrialMuscle(report.listRemove) = true;

    % Check if all trials are still very noisy
    if all(badTrialMuscle)
        warning('All trials are still full of EMG activity. Consider excluding this participant.'); return;
    elseif any(badTrialMuscle)
        EEG = pop_rejepoch(EEG,badTrialMuscle,0);
    else
        fprintf('Nice, there are no very contaminated trials for rejection.\n');
    end

else
    fprintf('Nice, no EMG found in the data. Skipping trial/channel interpolation.\n');
    badElecsPerTrial = cell(1,EEG.trials);
    report = [];
    report.listFixed  = [];
    report.listRemove = [];
end

% 2. Remove all
% BadEpochs = sum(slopesChannelsxEpochs, 1);
% badTrialMuscleTreshhold = 0;
% badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;
%
% while sum(badTrialMuscle) > 0.5*EEG.trials
%     warning('Increasing the minimum number of allowed EMG-contaminated channels: %d -> %d!',badTrialMuscleTreshhold,badTrialMuscleTreshhold+1);
%     badTrialMuscleTreshhold = badTrialMuscleTreshhold+1;
%     badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;
% end

% Log
EEG.ALSUTRECHT.epochRejections.InterpTrialInfo = badElecsPerTrial;
EEG.ALSUTRECHT.epochRejections.InterpReport    = report;
EEG.ALSUTRECHT.epochRejections.interpEpochs    = length(report.listFixed);

% Note the number of trials
NumberTrials(3) = EEG.trials;

% =============================================
% C. Detection using variance and the G-ESD method
fprintf('\n--------------------------------\n');
fprintf('G-ESD method with variance\n');
fprintf('--------------------------------\n');
% EEG = do_gsd(EEG);
fprintf('Not done.\n');

% Note the number of trials
NumberTrials(4) = EEG.trials;

% Log
EEG.ALSUTRECHT.epochRejections.initialEpochs       = NumberTrials(1);
EEG.ALSUTRECHT.epochRejections.afterEEGLABEpochs   = NumberTrials(2);
EEG.ALSUTRECHT.epochRejections.afterEMGSlopeEpochs = NumberTrials(3);
EEG.ALSUTRECHT.epochRejections.remainingEpochs     = NumberTrials(4);
EEG.ALSUTRECHT.epochRejections.proportionOfEpochsRejected = (EEG.ALSUTRECHT.epochRejections.initialEpochs - EEG.ALSUTRECHT.epochRejections.remainingEpochs) / EEG.ALSUTRECHT.epochRejections.initialEpochs;

end