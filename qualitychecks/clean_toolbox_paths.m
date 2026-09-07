function clean_toolbox_paths()
% CLEAN_TOOLBOX_PATHS Checks and removes FieldTrip, FieldTrip-lite, and EEGLAB
% paths from the active MATLAB search path.

% Split the current search path into individual folders
path_list = strsplit(path, pathsep);

% Target keywords to match (case-insensitive)
% 'fieldtrip' catches both standard FieldTrip and fieldtrip-lite
target_patterns = {'fieldtrip', 'eeglab'};

% Match folders containing any of the target patterns
match_idx = false(size(path_list));
for k = 1:numel(target_patterns)
    match_idx = match_idx | contains(path_list, target_patterns{k}, 'IgnoreCase', true);
end

target_paths = path_list(match_idx);

% Remove empty entries if present
target_paths = target_paths(~cellfun('isempty', target_paths));

% Remove from MATLAB path
if ~isempty(target_paths)
    % Join paths using the platform path separator to avoid argument limits
    rmpath(strjoin(target_paths, pathsep));
    fprintf('Removed %d path(s) related to FieldTrip / EEGLAB from the MATLAB path.\n', numel(target_paths));
else
    fprintf('No FieldTrip or EEGLAB paths found in the MATLAB path.\n');
end

end