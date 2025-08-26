function file_list = check_duplicates(root_path, folders_to_exclude)
%CHECK_DUPLICATES Recursively finds .m files and checks for duplicates.
%   FILE_LIST = CHECK_FOR_DUPLICATES(ROOT_PATH, FOLDERS_TO_EXCLUDE) searches
%   the folder ROOT_PATH for .m files, ignoring any files within subfolders
%   listed in FOLDERS_TO_EXCLUDE.
%
%   Arguments:
%   - root_path: Path to the directory you want to search.
%   - folders_to_exclude: (Optional) A cell array of strings with folder
%                         names to exclude, e.g., {'external', 'unused'}.
%
%   Returns:
%   - file_list: A structure array with details of the found files.

fprintf('\nChecking for duplicate functions...\n');

% --- 1. Validate Input ---
if nargin < 2
    folders_to_exclude = {}; % Default to an empty list if not provided
end
if ~isfolder(root_path)
    error('The provided path "%s" is not a valid folder.', root_path);
end

% --- 2. Find ALL .m Files Recursively ---
search_pattern = fullfile(root_path, '**', '*.m');
all_files = dir(search_pattern);

if isempty(all_files)
    warning('No .m files were found in the specified path.');
    file_list = [];
    return;
end

% --- 3. NEW: Filter out Excluded Folders ---
if ~isempty(folders_to_exclude)
    original_count = numel(all_files);

    % Create a logical index of which files to keep. Start with all true.
    keep_idx = true(1, original_count);

    for i = 1:original_count
        % Split the file's folder path by the separator ('\' or '/')
        path_parts = split(all_files(i).folder, filesep);

        % Check if any of the path parts is in our exclusion list.
        % The `intersect` function finds common elements between the two sets.
        if ~isempty(intersect(path_parts, folders_to_exclude))
            keep_idx(i) = false; % If a match is found, mark for removal.
        end
    end

    % Keep only the files that were not marked for removal.
    file_list = all_files(keep_idx);

    fprintf('Excluded %d files. Checking %d files...\n', ...
        original_count - numel(file_list), numel(file_list));
else
    file_list = all_files; % No filtering needed
    fprintf('Found %d .m files. Checking for duplicates...\n', numel(file_list));
end

if isempty(file_list)
    warning('No .m files remain after filtering.');
    return;
end

% --- 4. Check Each Remaining File for Duplicates ---
[file_list.path] = deal("");
[file_list.all_instances] = deal({});
[file_list.is_duplicate] = deal(false);

for i = 1:numel(file_list)
    file_list(i).path = fullfile(file_list(i).folder, file_list(i).name);
    instances_found = which(file_list(i).name, '-all');
    file_list(i).all_instances = instances_found;
    file_list(i).is_duplicate = (numel(instances_found) > 1);
end

% Check and report
if all(~[file_list.is_duplicate])
    fprintf('Duplicate check complete. All good!\n');
else
    error('Duplicate functions found!');
end

end