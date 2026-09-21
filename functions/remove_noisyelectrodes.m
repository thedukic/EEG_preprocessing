function [EEG, bad_channels_1] = remove_noisyelectrodes(EEG, cfg)
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

channel_mask = strcmp({EEG.chanlocs.type}, 'EEG');
channel_labels = {EEG.chanlocs(channel_mask).labels};
assert(all(channel_mask));

% =========================================================================
% PREP function
% =========================================================================
fprintf('Detecting noisy electrodes using the PREP toolbox...\n');

if cfg.channel.ransacOff
    disp('Deterministic method is used (ransac is off, thus iterations are not done).');

    % Run
    noisyOut = findNoisyChannels(EEG, cfg.channel);

    % Extract
    bad_channels_1 = noisyOut.noisyChannels.all;

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
    % bad_channels_1 = find(badness_percent >= cfg.channel.iter.frc);
end

% Report
num_removed_1 = length(bad_channels_1);
fprintf('PREP detected: %d\n', num_removed_1);

% =========================================================================
% Power spectra slope
% =========================================================================
% Determine how many we can still remove
channels_initial   = sum(strcmp({EEG.allchans.type}, 'EEG'));
channels_current   = sum(channel_mask) - num_removed_1;
max_remove_setting = round(cfg.channel.prop_badchan_max * channels_initial);
max_remove_here    = max_remove_setting - (channels_initial - channels_current);

fprintf('\nDetecting noisy electrodes using power slopes...\n');
if max_remove_here > 0
    % Estimate log-log slopes
    slopesChannelsxEpochs = detect_emg(EEG, cfg);

    % Detect noisy channels
    fprintf('Slope treshold: %1.2f\n', cfg.emg.slope_threshold_1);
    emgSlopeTimeAvg = mean(slopesChannelsxEpochs > cfg.emg.slope_threshold_1, 2);

    bad_channels_2    = find(emgSlopeTimeAvg > cfg.emg.slope_time);
    inital_number     = length(bad_channels_2);
    inital_proportion = inital_number / length(emgSlopeTimeAvg);

    if inital_number > max_remove_here
        fprintf('Warning: Too many electrodes (N = %d, max = %d) are marked for rejection based on their slope.\n', inital_number, max_remove_here);
        bad_channels_2sorted = sort(emgSlopeTimeAvg, 1, 'descend');
        emgSlopeTimeNew      = bad_channels_2sorted(max_remove_here, 1);
        bad_channels_2       = find(emgSlopeTimeAvg >= emgSlopeTimeNew);
        fprintf('Lowering that to N = %d.\n', length(bad_channels_2));
    end

    fprintf('Aberrant slope(s) detected: %d\n', length(bad_channels_2));

else
    warning('Skippping... Too many electrodes (N = %d, max = %d) have already been marked for rejection by the PREP toolbox.', num_removed_1, max_remove_setting);
    bad_channels_2   = [];
    inital_number     = NaN;
    inital_proportion = NaN;
end

% =========================================================================
% Combine abd remove (PREP and slope)
% =========================================================================
badElectrodes = unique([bad_channels_1(:); bad_channels_2(:)]);
fprintf('\nTotal detected: %d\n', length(badElectrodes));

% Remove
if ~isempty(badElectrodes)
    badElectrodes = channel_labels(badElectrodes);
    EEG = pop_select(EEG, 'nochannel', badElectrodes);
end

% Log
EEG.ALSUTRECHT.badchaninfo.prep.electrodes           = channel_labels(bad_channels_1);
EEG.ALSUTRECHT.badchaninfo.slope.electrodes          = channel_labels(bad_channels_2);
EEG.ALSUTRECHT.badchaninfo.slope.maxThatCanBeRemoved = max_remove_setting;
EEG.ALSUTRECHT.badchaninfo.slope.initalNumber        = inital_number;
EEG.ALSUTRECHT.badchaninfo.slope.initalProportion    = inital_proportion;

% =========================================================================
% Merge all bad electrodes
% =========================================================================
subfields = {'offsets', 'flat', 'prep', 'slope'};
bad_electrodes = [];

for i_field = subfields
    fn = i_field{1};
    if isfield(EEG.ALSUTRECHT.badchaninfo, fn) && ...
            isfield(EEG.ALSUTRECHT.badchaninfo.(fn), 'electrodes') && ...
            ~isempty(EEG.ALSUTRECHT.badchaninfo.(fn).electrodes)

        vals = EEG.ALSUTRECHT.badchaninfo.(fn).electrodes;
        if iscell(vals) && isempty(bad_electrodes)
            bad_electrodes = {};
        end
        bad_electrodes = [bad_electrodes, vals(:)']; %#ok<AGROW>
    end
end

if isempty(bad_electrodes)
    EEG.ALSUTRECHT.badchaninfo.badElectrodes = [];
else
    EEG.ALSUTRECHT.badchaninfo.badElectrodes = unique(bad_electrodes);
end

end