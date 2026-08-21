function [EEG, badElectrodes1] = remove_noisyelectrodes(EEG, cfg)
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

fprintf('\n================================\n');
fprintf('Detecting noisy electrodes\n');
fprintf('================================\n');

% % Detect which channels are EEG/EXT
% chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% chanext = {EEG.chanlocs(~chaneeg).labels};
% chanLabelsEEG = {EEG.chanlocs(chaneeg).labels};
%
% % Remove EXT channels, we dont want to check them
% EEGTMP = pop_select(EEG,'nochannel',chanext);

channel_mask = strcmp({EEG.chanlocs.type}, 'EEG');
channel_labels = {EEG.chanlocs(channel_mask).labels};
assert(all(channel_mask));

% =========================================================================
% PREP function
fprintf('Detecting noisy electrodes using the PREP toolbox...\n');

if cfg.channel.ransacOff
    disp('Deterministic method is used (ransac is off, thus iterations are not done).');

    % Run
    noisyOut = findNoisyChannels(EEG, cfg.channel);

    % Extract
    badElectrodes1 = noisyOut.noisyChannels.all;

else
    % badElectrodes_iter = false(EEG.nbchan,cfg.channel.iter.num);
    %
    % disp('Bad channel detection starting...');
    % for i = 1:cfg.channel.iter.num
    %     disp(['Iteration ' num2str(i) '/' num2str(cfg.channel.iter.num)]);
    %
    %     noisyOut = findNoisyChannels(EEG,cfg.channel);
    %     badElectrodes_iter(noisyOut.noisyChannels.all,i) = true;
    % end
    %
    % % RANSAC stability
    % badness_percent = sum(badElectrodes_iter,2) / size(badElectrodes_iter,2);
    %
    % % Check if there are too many bad channels detected (ie >iter.rejmax)
    % % If so, raise iter.frc to match iter.rejmax
    % if ~isempty(cfg.channel.iter.rejmax)
    %     if cfg.channel.iter.rejmax<1 && cfg.channel.iter.rejmax>0
    %         cfg.channel.iter.rejmax = round(EEG.nbchan * cfg.channel.iter.rejmax);
    %     end
    %     [sort_val, ~ ] = sort(badness_percent, 'descend');
    %
    %     if sort_val(cfg.channel.iter.rejmax) > cfg.channel.iter.frc
    %         cfg.channel.iter.frc = sort_val(cfg.channel.iter.rejmax);
    %     end
    % end
    %
    % % Final detection
    % badElectrodes = find(badness_percent >= cfg.channel.iter.frc);
end

% Report
Nremoved1 = length(badElectrodes1);
fprintf('PREP detected: %d\n', Nremoved1);

% =========================================================================
% Determine how many we can still remove
totalInitialChannels = sum(strcmp({EEG.allchans.type}, 'EEG'));
currentChannels      = sum(channel_mask) - Nremoved1;
maxThatCanBeRemoved  = round(cfg.channel.prop_badchan_max * totalInitialChannels);
youCanRejectThisManyChannelsHere = maxThatCanBeRemoved - (totalInitialChannels - currentChannels);

% =========================================================================
% Power spectra slope
fprintf('\nDetecting noisy electrodes using power slopes...\n');
if youCanRejectThisManyChannelsHere > 0
    % Estimate log-log slopes
    slopesChannelsxEpochs = detect_emg(EEG, cfg);

    % Detect noisy channels
    fprintf('Slope treshold: %1.2f\n', cfg.emg.slope_threshold_1);
    emgSlopeTimeAvg = mean(slopesChannelsxEpochs > cfg.emg.slope_threshold_1, 2);

    badElectrodes2   = find(emgSlopeTimeAvg > cfg.emg.slope_time);
    initalNumber     = length(badElectrodes2);
    initalProportion = initalNumber / length(emgSlopeTimeAvg);

    if initalNumber > youCanRejectThisManyChannelsHere
        fprintf('Warning: Too many electrodes (N = %d, max = %d) are marked for rejection based on their slope.\n', initalNumber, youCanRejectThisManyChannelsHere);
        badElectrodes2sorted = sort(emgSlopeTimeAvg, 1, 'descend');
        emgSlopeTimeNew      = badElectrodes2sorted(youCanRejectThisManyChannelsHere, 1);
        badElectrodes2       = find(emgSlopeTimeAvg >= emgSlopeTimeNew);
        fprintf('Lowering that to N = %d.\n', length(badElectrodes2));
    end

    fprintf('Aberrant slope(s) detected: %d\n', length(badElectrodes2));

else
    warning('Skippping... Too many electrodes (N = %d, max = %d) have already been marked for rejection by the PREP toolbox.', Nremoved1, maxThatCanBeRemoved);
    badElectrodes2   = [];
    initalNumber     = NaN;
    initalProportion = NaN;
end

% =========================================================================
% Combine
badElectrodes = unique([badElectrodes1(:); badElectrodes2(:)]);
fprintf('\nTotal detected: %d\n', length(badElectrodes));

% Remove
if ~isempty(badElectrodes)
    badElectrodes = channel_labels(badElectrodes);
    EEG = pop_select(EEG, 'nochannel', badElectrodes);
end

% Log
EEG.ALSUTRECHT.badchaninfo.prep.electrodes         = channel_labels(badElectrodes1);
EEG.ALSUTRECHT.badchaninfo.slope.electrodes          = channel_labels(badElectrodes2);
EEG.ALSUTRECHT.badchaninfo.slope.maxThatCanBeRemoved = maxThatCanBeRemoved;
EEG.ALSUTRECHT.badchaninfo.slope.initalNumber        = initalNumber;
EEG.ALSUTRECHT.badchaninfo.slope.initalProportion    = initalProportion;

end