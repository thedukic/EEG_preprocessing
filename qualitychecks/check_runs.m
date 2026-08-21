function list_failed = check_runs(myPaths)

% Initialise the list
list_failed = [];

% -------------------------------------------------------------------------
% Initialise a flag to track if any folders are missing
allProcessed = true;

% Get a list of all subfolders within the current path.
dirInfo = dir(myPaths.preproc);
subfolders = {dirInfo([dirInfo.isdir] & ~strcmp({dirInfo.name}, '.') & ~strcmp({dirInfo.name}, '..')).name};
subfolders = subfolders(contains(subfolders, 'ALS'));

num_folders = length(subfolders);
num_subjects = length(myPaths.subjects);
fprintf('%d participants should have been processed.\n', num_subjects);
fprintf('%d folders found.\n', num_folders);

% mask_diff = setdiff(myPaths.subjects, subfolders);

for j = 1:num_folders
    subfolderPath = fullfile(myPaths.preproc, subfolders{j}, '2', 'data');

    % Check for "_1b.mat" files in the current subfolder.
    fileList = dir(fullfile(subfolderPath, '*_b.mat'));

    if isempty(fileList)
        fprintf('%d. %s: processing failed.\n', j, subfolders{j});
        list_failed = [list_failed, j];
        allProcessed = false;
    end
end

% -------------------------------------------------------------------------
% Initialise a flag to track if any folders are missing
allPresent = true;

% Loop through each subject folder in the cell array
for i = 1:num_subjects
    % Extract the folder name from the cell
    folderName = myPaths.subjects{i};

    % Construct the full expected path
    fullPath = fullfile(myPaths.preproc, folderName);

    % Check if the folder exists
    if ~exist(fullPath, 'dir')
        fprintf('Missing folder: %s\n', fullPath);
        allPresent = false;
        list_failed = [list_failed i];
    end
end

% Make sure that there are not duplicates
list_failed = unique(list_failed);

% Report
if allProcessed & allPresent
    fprintf('All participants are successfully processed.\n');
else
    fprintf('Not all particiapnts are processed (N = %d).\n', length(list_failed));
end

end