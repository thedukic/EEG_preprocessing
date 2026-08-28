function DATA = fix_events(DATA)
% FIX_EVENTS Converts numeric event types to char/string.
% Implemented to harmonise Dublin EO data with the expected string format.

fprintf('\n================================\n');
fprintf('Fixing events (numeric -> char)\n');
fprintf('================================\n');

num_blocks = length(DATA);

for i_block = 1:num_blocks
    if isempty(DATA(i_block).event) || ~isfield(DATA(i_block).event, 'type')
        fprintf('Block %d: No events found.\n', i_block);
        continue;
    end

    % Check if event types are numeric
    if isnumeric(DATA(i_block).event(1).type)
        fprintf('Block %d:\n', i_block);

        % Extract all numeric types
        raw_types = [DATA(i_block).event.type];

        % Find unique numeric values for logging
        unique_types = unique(raw_types);

        % Convert all numeric types to strings in one vectorised step
        str_types = arrayfun(@num2str, raw_types, 'UniformOutput', false);

        % Reassign back to the event structure
        [DATA(i_block).event.type] = str_types{:};

        % Sync urevent if it exists
        if isfield(DATA(i_block), 'urevent') && ~isempty(DATA(i_block).urevent)
            raw_ur_types = [DATA(i_block).urevent.type];
            if isnumeric(raw_ur_types(1))
                str_ur_types = arrayfun(@num2str, raw_ur_types, 'UniformOutput', false);
                [DATA(i_block).urevent.type] = str_ur_types{:};
            end
        end

        % Log only the unique conversions
        for i_u = 1:length(unique_types)
            fprintf('Fixed: %d -> %s\n', unique_types(i_u), num2str(unique_types(i_u)));
        end
    else
        fprintf('Block %d: Not needed!\n', i_block);
    end
end

% Double-check consistency
DATA = eeg_checkset(DATA, 'eventconsistency');
DATA = eeg_checkset(DATA);

fprintf('\nDone!\n');
end