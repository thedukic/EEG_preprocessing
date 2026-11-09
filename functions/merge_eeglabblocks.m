function EEG = merge_eeglabblocks(EEG)
% MERGE_EEGLABBLOCKS
% Merges multiple continuous EEGLAB dataset blocks into a single continuous struct.
% Removes pre-existing edge boundary markers (latency <= 1 or >= pnts) prior to
% merging to prevent duplicate or out-of-bounds boundary events in pop_mergeset.

fprintf('\n================================\n');
fprintf('Merging blocks of data\n');
fprintf('================================\n');

num_block = length(EEG);
block_duration = NaN(num_block, 1);

for i_block = 1:num_block
    % Record individual block duration in seconds
    block_duration(i_block) = size(EEG(i_block).data, 2) / EEG(i_block).srate;

    if isfield(EEG(i_block), 'event') && ~isempty(EEG(i_block).event)
        n_pnts = EEG(i_block).pnts;

        % 1. Remove leading boundary event (latency <= 1.0 or == 0.5)
        if ~isempty(EEG(i_block).event) && ...
                (ischar(EEG(i_block).event(1).type) || isstring(EEG(i_block).event(1).type)) && ...
                strcmpi(EEG(i_block).event(1).type, 'boundary') && ...
                EEG(i_block).event(1).latency <= 1.0
            EEG(i_block).event(1) = [];
        end

        % 2. Remove trailing boundary event (latency >= pnts or >= pnts - 0.5)
        % Note: Re-checking ~isempty prevents a crash if the block only had 1
        % event that was just deleted in the step above.
        if ~isempty(EEG(i_block).event) && ...
                (ischar(EEG(i_block).event(end).type) || isstring(EEG(i_block).event(end).type)) && ...
                strcmpi(EEG(i_block).event(end).type, 'boundary') && ...
                EEG(i_block).event(end).latency >= (n_pnts - 0.5)
            EEG(i_block).event(end) = [];
        end
    end
end

% Handle single block case
if num_block == 1
    fprintf('Only 1 block provided. Merging skipped.\n');
    EEG.ALSUTRECHT.block.block_duration = block_duration;
    EEG = eeg_checkset(EEG);
    fprintf('Done!\n');
    return;
end

% Harmonise struct fields across all blocks to prevent pop_mergeset crash
% (pop_mergeset fails if one event struct has different fields than another)
all_fields = {};
for i_block = 1:num_block
    if isfield(EEG(i_block), 'event') && ~isempty(EEG(i_block).event)
        all_fields = union(all_fields, fieldnames(EEG(i_block).event), 'stable');
    end
end
if isempty(all_fields)
    all_fields = {'type', 'latency', 'duration'};
end

% Create a dummy struct containing all necessary fields
dummy_struct = cell2struct(cell(length(all_fields), 1), all_fields, 1);
dummy_struct.type     = 'dummy';
dummy_struct.latency  = 1;
dummy_struct.duration = 1;

% 1. Add a dummy event to every dataset to ensure consistency
for i_block = 1:num_block
    if isempty(EEG(i_block).event)
        EEG(i_block).event = dummy_struct;
    else
        % Ensure any missing fields are padded into the existing struct
        curr_fields = fieldnames(EEG(i_block).event);
        missing_fields = setdiff(all_fields, curr_fields);
        for f = 1:length(missing_fields)
            [EEG(i_block).event.(missing_fields{f})] = deal([]);
        end

        EEG(i_block).event(end+1) = dummy_struct;
    end
end

% 2. Merge the datasets
EEG = pop_mergeset(EEG, 1:num_block);

% 3. Identify and remove all 'dummy' events from the merged dataset
% (num2str prevents strcmp from crashing on numeric hardware triggers)
ev_types = cellfun(@num2str, {EEG.event.type}, 'UniformOutput', false);
dummy_idx = strcmpi(ev_types, 'dummy');
EEG.event(dummy_idx) = [];

% 4. Clear desynchronised urevent table to force clean regeneration
EEG.urevent = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% Validate that boundary events were inserted at concatenation points
ev_types = cellfun(@num2str, {EEG.event.type}, 'UniformOutput', false);
n_boundaries = sum(strcmpi(ev_types, 'boundary'));
assert(n_boundaries >= (num_block - 1), 'Expected at least %d boundary events after merging, found %d.', num_block - 1, n_boundaries);

% Store block metadata
EEG.ALSUTRECHT.block.block_duration = block_duration;

% Ensure internal struct consistency
EEG = eeg_checkset(EEG);

fprintf('Successfully merged %d blocks (%d total boundary events).\n', num_block, n_boundaries);
fprintf('Done!\n');

end