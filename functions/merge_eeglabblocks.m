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
                strcmpi(EEG(i_block).event(1).type, 'boundary') && ...
                EEG(i_block).event(1).latency <= 1.0
            EEG(i_block).event(1) = [];
        end

        % 2. Remove trailing boundary event (latency >= pnts or >= pnts - 0.5)
        if ~isempty(EEG(i_block).event) && ...
                strcmpi(EEG(i_block).event(end).type, 'boundary') && ...
                EEG(i_block).event(end).latency >= (n_pnts - 0.5)
            EEG(i_block).event(end) = [];
        end
    end
end

% Handle single block case
if num_block == 1
    fprintf('Only 1 block provided. Merging skipped.\n');
    EEG.ALSUTRECHT.blockinfo.block_duration = block_duration;
    EEG = eeg_checkset(EEG);
    fprintf('Done!\n');
    return;
end

% Merge datasets
EEG = pop_mergeset(EEG, 1:num_block);

% Validate that boundary events were inserted at concatenation points
if isfield(EEG, 'event') && ~isempty(EEG.event)
    event_types = {EEG.event.type};
    n_boundaries = sum(strcmpi(event_types, 'boundary'));
    assert(n_boundaries >= (num_block - 1), ...
        'Expected at least %d boundary events after merging, found %d.', ...
        num_block - 1, n_boundaries);
end

% Store block metadata
EEG.ALSUTRECHT.blockinfo.block_duration = block_duration;

% Ensure internal struct consistency
EEG = eeg_checkset(EEG);

fprintf('Successfully merged %d blocks (%d total boundary events).\n', num_block, n_boundaries);
fprintf('Done!\n');

end