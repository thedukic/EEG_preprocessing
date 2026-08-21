function [list_blocks, num_blocks] = list_datasets(datapath,thistask)

fprintf('\n================================\n');
fprintf('Looking up EEG files\n');
fprintf('================================\n');

if strcmpi(thistask, 'RS')
    % Resting-state
    expected_file_count = 3;

    % 1. Search for files and check the count for each task
    % The tmp variables now hold the dir structure AND the check is done.
    tmp1 = checkAndGetFiles(datapath, '*EO*', expected_file_count);
    tmp2 = checkAndGetFiles(datapath, '*EC*', expected_file_count);

    % 2. Combine the names
    % This is still done in the main script to consolidate all your data paths
    dataname = {tmp1.name, tmp2.name};

elseif strcmpi(thistask, 'MT')
    % Motor tasks
    expected_file_count = 1;

    % 1. Search for files and check the count for each task
    % The tmp variables now hold the dir structure AND the check is done.
    tmp1 = checkAndGetFiles(datapath, '*MT2*', expected_file_count);
    tmp2 = checkAndGetFiles(datapath, '*MT3*', expected_file_count);
    % tmp3 = checkAndGetFiles(datapath, '*MT5*', expected_file_count);

    % 2. Combine the names
    % This is still done in the main script to consolidate all your data paths
    % dataname = {tmp1.name, tmp2.name, tmp3.name};
    dataname = {tmp1.name, tmp2.name};

else
    % SART / MMN / (EO / EC)
    expected_file_count = 3;

    % The tmp variables now hold the dir structure AND the check is done.
    searchPattern = ['*' thistask '*'];
    tmp1 = checkAndGetFiles(datapath, searchPattern, expected_file_count);

    % 2. Combine the names
    % This is still done in the main script to consolidate all your data paths
    dataname = {tmp1.name};

end

list_blocks = fullfile(datapath,dataname);
num_blocks = length(list_blocks);

end

% =========================================================================
% Helper function
% =========================================================================
function tmpStruct = checkAndGetFiles(datapath, searchPattern, expectedCount)
% CHECKANDGETFILES searches for files, checks the count, and returns the dir structure.
%
%   TMPSTRUCT = CHECKANDGETFILES(DATAPATH, SEARCHPATTERN, EXPECTEDCOUNT)
%   Performs a 'dir' search, checks if the resulting file count matches
%   EXPECTEDCOUNT, issues a warning if they do not match, and returns the
%   resulting structure array.
%
%   Inputs:
%     datapath:      The full path to the directory (e.g., 'C:\Data\Subject1').
%     searchPattern: The pattern to search for (e.g., '*MT2*').
%     expectedCount: The exact number of files expected (e.g., 1).
%
%   Output:
%     tmpStruct:     The structure array returned by the 'dir' function.

% 1. Perform the search
fullSearchPattern = [datapath filesep searchPattern]; % Use filesep for cross-platform compatibility
tmpStruct = dir(fullSearchPattern);

% 2. Check the count
actualCount = numel(tmpStruct);

% Issue a warning with details
if actualCount ~= expectedCount
    warning('Expected %d file(s) for pattern "%s" in path "%s", but found %d.', ...
        expectedCount, searchPattern, datapath, actualCount);

    % Create a cell array of names and format for printing
    fileNamesCell = {tmpStruct.name};

    % Print the success message and the list in one go
    fprintf('Great, found %d file(s) for pattern "%s":\n', expectedCount, searchPattern);
    fprintf('- %s\n', fileNamesCell{:});
end

end