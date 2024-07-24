function [EEG, badElectrodes1] = remove_noisyelec(EEG,cfgbch)
%
% In this function, however, instead "findNoisyChannels" is called (from PREP pipeline)
% Note: "findNoisyChannels" must be edited so that it is undeterministic!
%
% PREP pipeline:
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
% http://vislab.github.io/EEG-Clean-Tools/
%
% Input EEG data should be:
% - EEGLAB struct
% - Highpass filtered >0.5 Hz
% - Lowpass filtering is okay if > 80 Hz
% - Referenced, probably the best to the common-average
%
% =========================================================================
% SDukic edit, Feb 2024
% =========================================================================

% =========================================================================
% 1. PREP function
fprintf('\nDetecting and removing noisy electrodes using PREP toolbox...\n');

% Remove external channels
chanext = {EEG.chanlocs(~strcmp({EEG.chanlocs.type},'EEG')).labels};
EEGTMP  = pop_select(EEG,'nochannel',chanext);

if cfgbch.ransacOff
    disp('Deterministic method is used (ransac is off), thus iterations are not done.');
    noisyOut = findNoisyChannels(EEGTMP,cfgbch);

    badElectrodes1 = noisyOut.noisyChannels.all;
else
    % badElectrodes_iter = false(EEG.nbchan,cfgbch.iter.num);
    %
    % disp('Bad channel detection starting...');
    % for i = 1:cfgbch.iter.num
    %     disp(['Iteration ' num2str(i) '/' num2str(cfgbch.iter.num)]);
    %
    %     noisyOut = findNoisyChannels(EEG,cfgbch);
    %     badElectrodes_iter(noisyOut.noisyChannels.all,i) = true;
    % end
    %
    % % RANSAC stability
    % badness_percent = sum(badElectrodes_iter,2) / size(badElectrodes_iter,2);
    %
    % % Check if there are too many bad channels detected (ie >iter.rejmax)
    % % If so, raise iter.frc to match iter.rejmax
    % if ~isempty(cfgbch.iter.rejmax)
    %     if cfgbch.iter.rejmax<1 && cfgbch.iter.rejmax>0
    %         cfgbch.iter.rejmax = round(EEG.nbchan * cfgbch.iter.rejmax);
    %     end
    %     [sort_val, ~ ] = sort(badness_percent, 'descend');
    %
    %     if sort_val(cfgbch.iter.rejmax) > cfgbch.iter.frc
    %         cfgbch.iter.frc = sort_val(cfgbch.iter.rejmax);
    %     end
    % end
    %
    % % Final detection
    % badElectrodes = find(badness_percent >= cfgbch.iter.frc);
end

% 2. Detect EMG-contaminated channels
% Estimate log-log power spectra
slopesChannelsxEpochs = dected_emg(EEGTMP);

% Detect noisy channels
badElectrodes2 = mean(slopesChannelsxEpochs>cfgbch.muscleSlopeThreshold,2);
badElectrodes2 = find(badElectrodes2>cfgbch.MuscleSlopeTime);

% Log
EEG.ALSUTRECHT.badchaninfo.PREPElectrodes = {EEG.chanlocs(badElectrodes1).labels};
EEG.ALSUTRECHT.badchaninfo.EMGSlope       = {EEG.chanlocs(badElectrodes2).labels};

% 3. Combine
badElectrodes = unique([badElectrodes1(:); badElectrodes2(:)]);

% 4. Remove
if ~isempty(badElectrodes)
    badElectrodes = {EEG.chanlocs(badElectrodes).labels};
    EEG = pop_select(EEG,'nochannel',badElectrodes);
end

% Report
fprintf('PREP detected %d bad electrodes.\n',length(badElectrodes1));
fprintf('EMG-slope detected %d bad electrodes.\n',length(badElectrodes2));

% ADD THIS
% RELAX_excluding_channels_and_epoching

end