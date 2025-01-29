function [EEG, badElectrodes1] = remove_noisyelectrodes(EEG,cfgbch)
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

% Detect which channels are EEG/EXT
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
chanext = {EEG.chanlocs(~chaneeg).labels};
chanLabelsEEG = {EEG.chanlocs(chaneeg).labels};

% Remove EXT channels, we dont want to check them
EEGTMP = pop_select(EEG,'nochannel',chanext);

% =========================================================================
% PREP function
fprintf('\nDetecting and removing noisy electrodes using PREP toolbox...\n');

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

% Report
Nremoved1 = length(badElectrodes1);
fprintf('PREP detected %d bad electrodes.\n',Nremoved1);

% =========================================================================
% Determine how many we can still remove
totalInitialChannels = sum(strcmp({EEG.allchans.type},'EEG'));
currentChannels      = sum(chaneeg) - Nremoved1;
maxThatCanBeRemoved  = round(cfgbch.maxProportionOfBadElec * totalInitialChannels);
youCanRejectThisManyChannelsHere = maxThatCanBeRemoved - (totalInitialChannels - currentChannels);

% =========================================================================
% Power spectra slope
if youCanRejectThisManyChannelsHere > 0
    % Estimate log-log slopes
    slopesChannelsxEpochs = detect_emg(EEGTMP,cfgbch);

    % Detect noisy channels
    muscleSlopeTimeAvg = mean(slopesChannelsxEpochs > cfgbch.muscleSlopeThreshold,2);

    badElectrodes2   = find(muscleSlopeTimeAvg > cfgbch.muscleSlopeTime);
    initalNumber     = length(badElectrodes2);
    initalProportion = initalNumber / length(muscleSlopeTimeAvg);

    if initalNumber > youCanRejectThisManyChannelsHere % initalProportion > cfgbch.maxProportionOfBadElec
        warning('Too many electrodes (N = %d, max = %d) are marked for rejection based on their EMG-slope.',initalNumber,youCanRejectThisManyChannelsHere);
        badElectrodes2sorted = sort(muscleSlopeTimeAvg,1,'descend');
        muscleSlopeTimeNew   = badElectrodes2sorted(youCanRejectThisManyChannelsHere,1);
        badElectrodes2       = find(muscleSlopeTimeAvg>=muscleSlopeTimeNew);
        warning('Lowering that to N = %d.',length(badElectrodes2));
    end

    fprintf('EMG-slope detected %d bad electrodes.\n',length(badElectrodes2));

else
    warning('Too many electrodes (N = %d, max = %d) are already marked for rejection by PREP toolbox.',Nremoved1,maxThatCanBeRemoved);
    badElectrodes2   = [];
    initalNumber     = NaN;
    initalProportion = NaN;
end

% =========================================================================
% Combine
badElectrodes = unique([badElectrodes1(:); badElectrodes2(:)]);

% Remove
if ~isempty(badElectrodes)
    badElectrodes = chanLabelsEEG(badElectrodes);
    EEG = pop_select(EEG,'nochannel',badElectrodes);
end

% Log
EEG.ALSUTRECHT.badchaninfo.maxThatCanBeRemoved      = maxThatCanBeRemoved;
EEG.ALSUTRECHT.badchaninfo.PREPElectrodes           = chanLabelsEEG(badElectrodes1);
EEG.ALSUTRECHT.badchaninfo.EMGSlope                 = chanLabelsEEG(badElectrodes2);
EEG.ALSUTRECHT.badchaninfo.EMGSlopeInitalSum        = initalNumber;
EEG.ALSUTRECHT.badchaninfo.EMGSlopeInitalProportion = initalProportion;

end