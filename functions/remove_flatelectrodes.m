function [EEG, badElectrodes] = remove_flatelectrodes(EEG,cfg)
%
% Removes flat electrodes (electrodes that were taken out)
% The code is based on the function from EEGLAB
% EEG = clean_flatlines(EEG,cfgbch.flatdur);
% But this code is data-driven as it estimates the recorded activity in each electrode
% And then, it estimates what could be an appropirate treshold for flat electrode detection
%
% Checked using:
% ALS37840 ALS T1: B2
% ALS37793 CON T1: A32

fprintf('\nChecking if there are any flat electrodes...\n');

% Minimum "flat" duration in [s]
if isempty(cfg.flatDuration)
    cfg.flatDuration = 4;
end

% Minimum "flat" duration in [samples]
T0 = cfg.flatDuration*EEG(1).srate;

% EXT/EMG electrodes are not checked (they should never be flat)
% otherchan = {EEG(1).chanlocs(~strcmp({EEG(1).chanlocs.type},'EEG')).labels};
% eegchan   = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EEG')).labels};
% EEGTMP    = pop_select(EEG,'nochannel',otherchan);
eegchan       = strcmp({EEG(1).chanlocs.type},'EEG');
eegchanLabels = {EEG(1).chanlocs(eegchan).labels};

NCHN = sum(eegchan);
NBLK = length(EEG);

% Remove possible slow drifts and strong 50 Hz noise
% These are likley to mask if the electrode is truly flat
% Especially the 50 Hz noise present in hospitals
fprintf('Temporarily filtering the EEG electrodes for better detection of flat electrodes.\n');
EEGTMP = filter_signal(EEG,[20 8],[1 8],1:NCHN,'eeglab');

fprintf('\n');
badElectrodes = cell(NBLK,1);
for i = 1:NBLK
    % This works better than the EEGLAB's manual treshold
    dataDiff = abs(diff(EEGTMP(i).data(1:NCHN,:),1,2));
    dataDiffMedian = median(dataDiff,2);
    % DiffMedianTreshold = min(dataDiffMedian);
    % DiffMedianTreshold = median(dataDiffMedian)
    DiffMedianTreshold = prctile(dataDiffMedian,5);
    fprintf('Block %d: The estimated treshold is %1.2f.\n', i,DiffMedianTreshold);

    % EEGLAB
    badElectrodesTmp = false(1,NCHN);
    for j = 1:NCHN
        zero_intervals = reshape(find(diff([false dataDiff(j,:)<=DiffMedianTreshold false])),2,[])';
        if max(zero_intervals(:,2) - zero_intervals(:,1)) > T0
            badElectrodesTmp(j) = true;
        end
    end

    % Report
    if any(badElectrodesTmp)
        badElectrodes{i} = eegchanLabels(badElectrodesTmp);
        fprintf('Block %d: Flat electrodes are found (N = %d)!\n', i,length(badElectrodes{i}));

        badElectrodesTmp = find(badElectrodesTmp);
        for j = 1:length(badElectrodes{i})
            fprintf('Electrode %s: Median voltage change is %1.2f.\n', badElectrodes{i}{j},dataDiffMedian(badElectrodesTmp(j)));
        end
        fprintf('\n');
    else
        fprintf('Block %d: Flat electrodes are not found.\n', i);
    end
end

% Check if all blocks have the same number of detected flat electrodes
N = cell2mat(cellfun(@(x) length(x),badElectrodes,'UniformOutput',false));
if range(N)~=0
    warning('Strange, not all blocks have the same number of detected flat electrodes.');
    flagDiscrep = true;
else
    flagDiscrep = false;
end

% Remove all flat electrodes
badElectrodes = unique(cat(2,badElectrodes{:}));
if ~isempty(badElectrodes)
    EEG = pop_select(EEG,'nochannel',badElectrodes);
    fprintf('The flat electrodes (N = %d) are removed from all data blocks (N = %d)!\n',length(badElectrodes),NBLK);
else
    badElectrodes = {};
    fprintf('Flat electrodes are not found in the dataset.\n');
end

% Log info
for i = 1:NBLK
    EEG(i).ALSUTRECHT.badchaninfo.flatElectrodes = badElectrodes;
    if flagDiscrep
        EEG(i).ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy = 1;
    else
        EEG(i).ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy = 0;
    end
end

% =========================================================================
% Detect gel birdges ???
% This might be inflated when lowpass filtering with low cutoff
% [bridge_pairs, EEG] = bemobil_find_gel_bridges(EEG,cfgbch.corrmax);
% bridge_pairs = NaN;

end