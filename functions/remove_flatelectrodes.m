function [EEG, elec_flat_labels] = remove_flatelectrodes(EEG, cfg)
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

fprintf('\n================================\n');
fprintf('Detecting flat electrodes\n');
fprintf('================================\n');

% Minimum "flat" duration in [samples]
T0 = cfg.channel.flatDuration * EEG(1).srate;

% EXT/EMG electrodes are not checked (they should never be flat)
eegchan       = strcmp({EEG(1).chanlocs.type}, 'EEG');
eegchanLabels = {EEG(1).chanlocs(eegchan).labels};

num_channel = sum(eegchan);
num_block = length(EEG);
fprintf('Checking %d channels across %d blocks.\n', num_channel, num_block);

% Remove possible slow drifts and strong 50 Hz noise
% These are likley to mask if the electrode is truly flat
% Especially the 50 Hz noise present in hospitals
fprintf('Temporarily filtering the EEG electrodes for better detection of flat electrodes.\n');
EEG_TMP = filter_signal(EEG, [20 8], [1 8], 1:num_channel, 'eeglab');

% Typically
elec_high_offset_indx = EEG(1).ALSUTRECHT.badchaninfo.offsets.electrodes_indx;
elec_high_offset = false(128, 1);
if ~isempty(elec_high_offset_indx)
    elec_high_offset_indx_tmp = elec_high_offset_indx;
    elec_high_offset_indx_tmp = elec_high_offset_indx_tmp(elec_high_offset_indx_tmp <= 128);
    elec_high_offset(elec_high_offset_indx_tmp) = true;
    fprintf('Electrodes with high offset were found earlier: %d\n', sum(elec_high_offset));
else
    fprintf('Electrodes with high offset were not found earlier\n');
end

fprintf('\n');
elec_flat_indx = cell(num_block,1);
elec_flat_labels = cell(num_block,1);

for i_block = 1:num_block
    data_diff = abs(diff(EEG_TMP(i_block).data(1:num_channel,:), 1, 2));
    data_diff_median = median(data_diff, 2);

    % The 5th percentile will flag CMS neighbours, but the '&' logic below protects them
    data_diff_median_threshold = prctile(data_diff_median, 5);

    fprintf('Block %d: The estimated threshold is %1.2f.\n', i_block, data_diff_median_threshold);

    % Detect low variance (flat) in the filtered data
    elec_flat_tmp = false(num_channel, 1);
    for i_channel = 1:num_channel
        zero_intervals = reshape(find(diff([false data_diff(i_channel,:) <= data_diff_median_threshold false])), 2, [])';
        if ~isempty(zero_intervals)
            if max(zero_intervals(:,2) - zero_intervals(:,1)) > T0
                elec_flat_tmp(i_channel) = true;
            end
        end
    end

    % Combine: A dangling electrode must have a massive offset AND a filtered flatline
    elec_flat_tmp = elec_high_offset(:) & elec_flat_tmp(:);

    % Report
    if any(elec_flat_tmp)
        elec_flat_labels{i_block} = eegchanLabels(elec_flat_tmp);
        fprintf('Block %d: Flat electrodes are found (N = %d)!\n', i_block, length(elec_flat_labels{i_block}));

        elec_flat_idx = find(elec_flat_tmp);
        for i_channel = 1:length(elec_flat_labels{i_block})
            ch_idx = elec_flat_idx(i_channel);
            fprintf('Electrode %s: High offset + Median voltage change %1.2f.\n', elec_flat_labels{i_block}{i_channel}, data_diff_median(ch_idx));
        end
        fprintf('\n');

        elec_flat_indx{i_block} = elec_flat_idx;
    else
        fprintf('Block %d: Flat electrodes are not found.\n', i_block);
    end
end

% Check if all blocks have the same number of detected flat electrodes
N = cell2mat(cellfun(@(x) length(x), elec_flat_labels, 'UniformOutput', false));
if range(N) ~= 0
    warning('Strange, not all blocks have the same number of detected flat electrodes.');
    flag_discrep = true;
else
    flag_discrep = false;
end

% Remove all flat electrodes
elec_flat_indx   = unique(cat(2, elec_flat_indx{:}));
elec_flat_labels = unique(cat(2, elec_flat_labels{:}));
if ~isempty(elec_flat_labels)
    % EEG = pop_select(EEG, 'nochannel', badElectrodes);
    % fprintf('The flat electrodes (N = %d) are removed from all data blocks (N = %d).\n', length(elec_flat_1), num_block);

    for i_block = 1:num_block
        EEG(i_block) = pop_interp(EEG(i_block), elec_flat_indx, 'spherical'); % !!!!!!!!!!!!!!!!!!!!!!!!!
    end
    fprintf('The flat electrodes (N = %d) are interpoalted in all data blocks (N = %d).\n', length(elec_flat_labels), num_block);
else
    elec_flat_labels = {};
    fprintf('Flat electrodes are not found in this dataset.\n');
end

% Log info
for i_block = 1:num_block
    EEG(i_block).ALSUTRECHT.badchaninfo.flat.electrodes = elec_flat_labels;
    if flag_discrep
        EEG(i_block).ALSUTRECHT.badchaninfo.flat.discrepancy = true;
    else
        EEG(i_block).ALSUTRECHT.badchaninfo.flat.discrepancy = false;
    end
end

end