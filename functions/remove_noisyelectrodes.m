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

% % Temporary filtering
% EEG.data = do_ref(EEG.data, EEG.nbchan);
% cut_off = 0.5;
% EEG_CHANNEL = do_filtering_method(EEG, cut_off, cfg);
EEG_CHANNEL = EEG;

% =========================================================================
% PREP function
% =========================================================================
fprintf('Detecting noisy electrodes using the PREP toolbox...\n');

% Run
noisyOut = findNoisyChannels(EEG_CHANNEL, cfg.channel.prep);

% Extract
bad_channels_1 = noisyOut.noisyChannels.all;

% Report
num_removed_1 = length(bad_channels_1);
fprintf('PREP detected: %d\n', num_removed_1);

% =========================================================================
% NEAR function
% =========================================================================
fprintf('Detecting noisy electrodes using the PREP toolbox...\n');

% Run
[EEG_TMP, flat_ch, lof_ch, periodo_ch, LOF_vec, thresh_lof_update] = NEAR_getBadChannels(EEG_CHANNEL, ...
    cfg.channel.near.isFlat, cfg.channel.near.flatWin, ...
    cfg.channel.near.isLOF, cfg.channel.near.thresh_lof, cfg.channel.near.dist_metric, cfg.channel.near.isAdapt, ...]
    cfg.channel.near.isPeriodogram, cfg.channel.near.frange, cfg.channel.near.winsize, cfg.channel.near.winov, cfg.channel.near.pthresh, 0);

% Extract
bad_channels_2 = sort(unique(union(union(flat_ch, lof_ch), periodo_ch)));

% Report
num_removed_2 = length(bad_channels_2);
fprintf('NEAR detected: %d\n', num_removed_2);

% =========================================================================
% Power spectra slope
% =========================================================================
% % Determine how many we can still remove
% channels_initial   = sum(strcmp({EEG.allchans.type}, 'EEG'));
% channels_current   = sum(channel_mask) - num_removed_1;
% max_remove_setting = round(cfg.channel.prop_badchan_max * channels_initial);
% max_remove_here    = max_remove_setting - (channels_initial - channels_current);
%
% fprintf('\nDetecting noisy electrodes using power slopes...\n');
% if max_remove_here > 0
%     % Estimate log-log slopes
%     slopesChannelsxEpochs = detect_emg(EEG, cfg);
%
%     % Detect noisy channels
%     fprintf('Slope treshold: %1.2f\n', cfg.emg.slope_threshold_1);
%     emgSlopeTimeAvg = mean(slopesChannelsxEpochs > cfg.emg.slope_threshold_1, 2);
%
%     bad_channels_2    = find(emgSlopeTimeAvg > cfg.emg.slope_time);
%     inital_number     = length(bad_channels_2);
%     inital_proportion = inital_number / length(emgSlopeTimeAvg);
%
%     if inital_number > max_remove_here
%         fprintf('Warning: Too many electrodes (N = %d, max = %d) are marked for rejection based on their slope.\n', inital_number, max_remove_here);
%         bad_channels_2sorted = sort(emgSlopeTimeAvg, 1, 'descend');
%         emgSlopeTimeNew      = bad_channels_2sorted(max_remove_here, 1);
%         bad_channels_2       = find(emgSlopeTimeAvg >= emgSlopeTimeNew);
%         fprintf('Lowering that to N = %d.\n', length(bad_channels_2));
%     end
%
%     fprintf('Aberrant slope(s) detected: %d\n', length(bad_channels_2));
%
% else
%     warning('Skippping... Too many electrodes (N = %d, max = %d) have already been marked for rejection by the PREP toolbox.', num_removed_1, max_remove_setting);
%     bad_channels_2   = [];
%     inital_number     = NaN;
%     inital_proportion = NaN;
% end

% =========================================================================
% Combine abd remove (PREP and slope)
% =========================================================================
clearvars EEG_CHANNEL EEG_TMP

electrode_bad = unique([bad_channels_1(:); bad_channels_2(:)]);
fprintf('\nTotal detected: %d\n', length(electrode_bad));

% Remove
if ~isempty(electrode_bad)
    electrode_bad = channel_labels(electrode_bad);
    EEG = pop_select(EEG, 'nochannel', electrode_bad);
end

% Log
EEG.ALSUTRECHT.bad_channel.prep.electrodes        = channel_labels(bad_channels_1);
EEG.ALSUTRECHT.bad_channel.prep.noisyOut          = noisyOut;
EEG.ALSUTRECHT.bad_channel.near.electrodes        = channel_labels(bad_channels_2);
EEG.ALSUTRECHT.bad_channel.near.LOF_vec           = LOF_vec;
EEG.ALSUTRECHT.bad_channel.near.thresh_lof_update = thresh_lof_update;
% EEG.ALSUTRECHT.bad_channel.slope.electrodes          = channel_labels(bad_channels_2);
% EEG.ALSUTRECHT.bad_channel.slope.maxThatCanBeRemoved = max_remove_setting;
% EEG.ALSUTRECHT.bad_channel.slope.initalNumber        = inital_number;
% EEG.ALSUTRECHT.bad_channel.slope.initalProportion    = inital_proportion;

% =========================================================================
% Merge all bad electrodes
% =========================================================================
subfields = {'offsets', 'flat', 'prep', 'near'};
bad_electrodes = [];

for i_field = subfields
    fn = i_field{1};
    if isfield(EEG.ALSUTRECHT.bad_channel, fn) && ...
            isfield(EEG.ALSUTRECHT.bad_channel.(fn), 'electrodes') && ...
            ~isempty(EEG.ALSUTRECHT.bad_channel.(fn).electrodes)

        vals = EEG.ALSUTRECHT.bad_channel.(fn).electrodes;
        if iscell(vals) && isempty(bad_electrodes)
            bad_electrodes = {};
        end
        bad_electrodes = [bad_electrodes, vals(:)']; %#ok<AGROW>
    end
end

if isempty(bad_electrodes)
    EEG.ALSUTRECHT.bad_channel.electrodes = [];
else
    EEG.ALSUTRECHT.bad_channel.electrodes = unique(bad_electrodes);
end
end

% function data = do_ref(data, num_channel)
% data_mean = sum(data, 1) / num_channel;
% data = data - data_mean;
% end
