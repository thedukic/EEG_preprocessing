function DATA = load_biosemidata(subject, cfg)
fprintf('\n================================\n');
fprintf('Loading data\n');
fprintf('================================\n');

% Load data
if strcmpi(subject.task, 'MT')
    % Motor tasks
    DATA = pop_biosig(subject.datablocks, 'channels', 1:168);
    assert(DATA(1).nbchan == 168, 'Unexpected number of channels');

    % fprintf('\nRemoving empty channels:\n');
    % Remove the unused EMG channels
    % 1. The two unused EMG channel for ADM (131-132)
    % 2. The rest of the unused EMG channels (143-160)
    % Loop through the datasets to avoid memory thrashing/throttling
    % for i = 1:length(DATA)
    %     DATA(i) = pop_select(DATA(i), 'rmchannel', [131:132, 143:160]);
    % end
    % % CMC tests
    % [cmc_topo, f_axis, cmc_all, cl95, fh] = compute_cmc_topography(DATA(1));
    % res = plot_cmc_rect_comparison(DATA(1));

    % Remove all EMG data as they will be processed separtely
    for i = 1:length(DATA)
        DATA(i) = pop_select(DATA(i), 'rmchannel', 129:160);
    end

else
    % Load EEG data
    DATA = pop_biosig(subject.datablocks);
    num_blocks = length(DATA);
    % Check if the number channels is as expected
    fprintf('\nChecking for empty electrodes in each block...\n');
    for i_block = 1:num_blocks
        switch DATA(i_block).nbchan
            case 136
                % As expected
            case 168
                % Utrect RS datasets recorded with Motor task settings
                % EEG = pop_biosig(subject.datablocks,'channels', 1:168);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 161:168]);
            case 264
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:264);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 257:264]);
            case 271
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:271);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 257:264]);
            case 143
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:143);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', 1:136);
            otherwise
                error('Unexpected number of channels!');
        end
    end
    assert(all([DATA(:).nbchan] == 136));
    fprintf('Done!\n');
end

% Ensure that data is with fs = 512 Hz
% Consider doing  the resampling (to 256 Hz) here
fprintf('\nChecking the sampling rate in each block...\n');
srate_all = [DATA(:).srate] > 512;
idx_resample = find(srate_all);
if ~isempty(idx_resample)
    % % Temporarily disable EEGLAB parallel processing
    % pop_editoptions('option_parallel', 0);

    % Loop through only the datasets that require downsampling
    for i = 1:length(idx_resample)
        idx = idx_resample(i);
        fprintf('Downsampling dataset %d of %d (from %d Hz to 512 Hz)...\n', i, length(idx_resample), DATA(idx).srate);
        DATA(idx) = pop_resample(DATA(idx), 512);
    end

    % % Re-enable parallel processing for subsequent steps
    % pop_editoptions('option_parallel', 1);
end
% Assert
assert(all([DATA(:).srate] == 512));
fprintf('Done!\n');


% Cut intruder blocks from other tasks
fprintf('\nChecking for task crosstalk (1 block = 2 tasks)...\n');
[DATA, subject.datablocks_trimmed] = trim_crosstask_data(DATA, subject, cfg);

end


function [DATA, flags_trimmed] = trim_crosstask_data(DATA, subject, cfg)
% TRIM_CROSSTASK_DATA Crops combined recording files containing multiple tasks.
%
% Syntax:
%   [DATA, flags_trimmed] = trim_crosstask_data(DATA, subject, cfg)
%
% Logic:
%   - If target events are at the START: crop from sample 1 to the first 255 after target events.
%   - If target events are at the END: crop from the last 255 before target events to the end of the file.

num_blocks = length(DATA);
flags_trimmed = false(1, num_blocks);

if ~isfield(cfg, 'trg') || ~isstruct(cfg.trg)
    return;
end

curr_task = lower(regexprep(subject.task, '\d+$', ''));
trg_fields = fieldnames(cfg.trg);

% -------------------------------------------------------------------------
% 1. Collect target and intruder trigger lists
% -------------------------------------------------------------------------
target_trgs       = [];
all_intruder_trgs = [];

for i_f = 1:length(trg_fields)
    fn = trg_fields{i_f};
    base_task = lower(regexprep(fn, '\d+$', ''));

    if startsWith(base_task, 'rs')
        % fprintf('No need! Resting-state data...'); continue;
        trg_vals = 999999;
    else
        trg_vals = cfg.trg.(fn){1};
    end

    if strcmp(base_task, curr_task)
        target_trgs = unique([target_trgs, trg_vals]);
    else
        all_intruder_trgs = unique([all_intruder_trgs, trg_vals]);
    end
end

if isempty(all_intruder_trgs) || isempty(target_trgs)
    fprintf('[Cross-Task Trim] %s: Good! Not needed.\n', subject.task);
    return;
end

% -------------------------------------------------------------------------
% 2. Inspect and crop each block
% -------------------------------------------------------------------------
margin_sec_255 = 0.5; % 0.5 s buffer relative to trigger 255

for i_block = 1:num_blocks
    block_str = sprintf('%s%d', subject.task, i_block);

    if isempty(DATA(i_block).event) || ~isfield(DATA(i_block).event, 'type')
        fprintf('[Cross-Task Trim] %s: Good! Not needed.\n', block_str);
        continue;
    end

    % Extract event types and latencies
    raw_types = {DATA(i_block).event.type};
    if isnumeric(raw_types{1})
        ev_types = [raw_types{:}];
    else
        ev_types = cellfun(@str2double, raw_types);
    end
    ev_latencies = [DATA(i_block).event.latency];

    valid_mask   = ~isnan(ev_types);
    ev_types     = ev_types(valid_mask);
    ev_latencies = ev_latencies(valid_mask);

    % Find target and intruder indices
    target_idx   = find(ismember(ev_types, target_trgs));
    intruder_idx = find(ismember(ev_types, all_intruder_trgs));

    % Only crop if both target and intruder triggers exist in this file
    if isempty(target_idx) || isempty(intruder_idx)
        fprintf('[Cross-Task Trim] %s: Good! Not needed.\n', block_str);
        continue;
    end

    fs = DATA(i_block).srate;
    pnts_total = DATA(i_block).pnts;
    margin_samples = round(margin_sec_255 * fs);

    min_target_lat   = min(ev_latencies(target_idx));
    max_target_lat   = max(ev_latencies(target_idx));
    min_intruder_lat = min(ev_latencies(intruder_idx));
    max_intruder_lat = max(ev_latencies(intruder_idx));

    % Case 1: Target task is strictly at the START (finished before intruder started)
    if max_target_lat < min_intruder_lat
        % Find first 255 after the target events
        % Returns only the first matching index
        idx_255 = find((ev_types == 255) & (ev_latencies > max_target_lat), 1, 'first');

        % if isempty(idx_255)
        %     idx_255 = find((ev_types == 255) & (ev_latencies >= min_intruder_lat));
        % end

        if ~isempty(idx_255)
            lat_255 = ev_latencies(idx_255);
            cut_end = lat_255 - margin_samples;
        else
            % cut_point = min_intruder_lat - margin_samples;
            error('check this');
        end

        cut_start = 1;
        pos_str   = 'START';

        % Case 2: Target task is strictly at the END (started after intruder finished)
    elseif min_target_lat > max_intruder_lat
        % Find the last 255 before target events begin
        idx_255 = find((ev_types == 255) & (ev_latencies < min_target_lat), 1, 'last');

        if ~isempty(idx_255)
            lat_255 = ev_latencies(idx_255);
            cut_start = lat_255 + margin_samples;
        else
            % % Fallback: cut 3.0 s after the intruder finished
            % cut_start = min(pnts_total, max_intruder_lat + round(margin_sec_backup * fs));
            error('check this');
        end

        cut_end = pnts_total;
        pos_str = 'END';

    else
        error('[Cross-Task Trim] %s: Triggers overlap or are interleaved. Skipping automated cut.', block_str);
    end

    % Apply crop
    orig_dur_s = pnts_total / fs / 60;
    new_dur_s  = (cut_end - cut_start + 1) / fs / 60;

    fprintf('[Cross-Task Trim] %s:\n', block_str);
    fprintf('   * Target task located at the %s of recording.\n', pos_str);
    fprintf('   * Cropping dataset: retaining samples %d to %d (%.1f min -> %.2f min)\n', ...
        cut_start, cut_end, orig_dur_s, new_dur_s);

    DATA(i_block) = pop_select(DATA(i_block), 'point', [cut_start, cut_end]);
    flags_trimmed(i_block) = true;

    % ---------------------------------------------------------------------
    % 3. Verify that no intruder triggers remain
    % ---------------------------------------------------------------------
    if ~isempty(DATA(i_block).event) && isfield(DATA(i_block).event, 'type')
        post_raw_types = {DATA(i_block).event.type};
        if isnumeric(post_raw_types{1})
            post_ev_types = [post_raw_types{:}];
        else
            post_ev_types = cellfun(@str2double, post_raw_types);
        end
        post_ev_types = post_ev_types(~isnan(post_ev_types));

        remaining_intruders = intersect(post_ev_types, all_intruder_trgs);
        assert(isempty(remaining_intruders), ...
            '[Cross-Task Trim] Block %d (%s) still contains intruder triggers after trimming: %s', ...
            i_block, block_str, mat2str(remaining_intruders));
    end
end

end