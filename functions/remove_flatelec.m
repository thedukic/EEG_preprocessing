function [EEG, badElectrodes] = remove_flatelec(EEG,cfg)
%
% Removes flat electrodes (electrodes that were taken out)
% The code is based on the function from EEGLAB
% EEG = clean_flatlines(EEG,cfgbch.flatdur);
% But this code is data-driven as it estimates the recorded activity in each channel
% And then, it estimates what could be an appropirate treshold for flat channel detection
%
% Checked using:
% ALS37840 ALS T1: B2
% ALS37793 CON T1: A32

fprintf('\nChecking if there are any flat channels...\n');

% Minimum "flat" duration in [s]
if isempty(cfg.flatDuration)
    cfg.flatDuration = 4;
end

% Minimum "flat" duration in [samples]
T0 = cfg.flatDuration*EEG(1).srate;

% EXT/EMG channels are not checked (they should never be flat)
otherchan = {EEG(1).chanlocs(~strcmp({EEG(1).chanlocs.type},'EEG')).labels};
eegchan   = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EEG')).labels};
EEGTMP    = pop_select(EEG,'nochannel',otherchan);

% Remove possible slow drifts and strong 50 Hz noise
% These are likley to mask if the channel is truly flat
% Especially the 50 Hz noise present in hospitals
fprintf('Temporary filtering for better detection of flat channels:\n');
EEGTMP = filter_signal(EEGTMP,[20 8],[1 8],1:length(eegchan),'eeglab');

NBLK = length(EEG);
NCHN = EEGTMP(1).nbchan;

badElectrodes = cell(NBLK,1);
for i = 1:NBLK
    % K = min(median(abs(diff(EEGTMP(i).data,1,2)),2));
    % K = median(median(abs(diff(EEGTMP(i).data,1,2)),2))
    K = prctile(median(abs(diff(EEGTMP(i).data,1,2)),2),10);

    % EEGLAB snippet
    removed_channels = false(1,NCHN);
    for j = 1:NCHN
        zero_intervals = reshape(find(diff([false abs(diff(EEGTMP(i).data(j,:)))<=K false])),2,[])';
        if max(zero_intervals(:,2) - zero_intervals(:,1)) > T0
            removed_channels(j) = true;
        end
    end

    % % My code: SDukic, August 2023
    % pow = NaN(NCHN,1);
    % T = floor(size(EEGTMP(i).data,2)/T0)*T0;
    % for j = 1:NCHN
    %     % pow(j) = median(var(reshape(EEGTMP(i).data(j,1:T),T0,[])));
    %     pow(j) = median(median(reshape(abs(EEGTMP(i).data(j,1:T)),T0,[])));
    % end
    % P = pow([1:64 97:NCHN]);
    % Z = (pow - mean(P)) / std(P);

    if any(removed_channels)
        badElectrodes{i} = eegchan(removed_channels);
        fprintf('BLOCK %d: Flat channels are found (N = %d)!\n', i, length(badElectrodes{i}));
    else
        fprintf('BLOCK %d: Flat channels are not found.\n', i);
    end
end

% Check if all blocks have the same number of detected flat channels
N = cell2mat(cellfun(@(x) length(x),badElectrodes,'UniformOutput',false));
if range(N)~=0
    warning('Strange, not all blocks have the same number of flat channels detected...');
    discrepflag = true;
else
    discrepflag = false;
end

% Remove all flat channels
badElectrodes = unique(cat(2,badElectrodes{:}));
if ~isempty(badElectrodes)
    EEG = pop_select(EEG,'nochannel',badElectrodes);
    fprintf('The flat channels (N = %d) are removed from all data blocks (N = %d)!\n',length(badElectrodes),NBLK);
    % disp(badElectrodes);
else
    badElectrodes = {};
end

% Log info
for i = 1:NBLK
    EEG(i).ALSUTRECHT.badchaninfo.flatElectrodes = badElectrodes;
    if discrepflag
        EEG(i).ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy = 1;
    else
        EEG(i).ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy = 0;
    end
end

% =========================================================================
% Gel birdges???
% This might be inflated when lowpass filtering with low cutoff
% [bridge_pairs, EEG] = bemobil_find_gel_bridges(EEG,cfgbch.corrmax);
% bridge_pairs = NaN;

end