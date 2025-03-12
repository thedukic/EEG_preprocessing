function check_failedruns(folderPaths)
% findFoldersWithout1bMat(folderPaths)
%
% Checks all folders within the given paths and identifies folders that do
% not contain a file ending with "_1b.mat".
%
% Inputs:
%   folderPaths: A cell array of strings, where each string is a folder path.
%
% Outputs:
%   Prints the paths of folders that lack a "_1b.mat" file.

if ~iscell(folderPaths)
    folderPaths = {folderPaths}; % Convert to cell if single string
end

for i = 1:length(folderPaths)
    currentPath = folderPaths{i};

    if ~isfolder(currentPath)
        fprintf('Warning: Path "%s" is not a folder.\n', currentPath); continue;
    end

    % Get a list of all subfolders within the current path.
    dirInfo = dir(currentPath);
    subfolders = {dirInfo([dirInfo.isdir] & ~strcmp({dirInfo.name}, '.') & ~strcmp({dirInfo.name}, '..')).name};

    for j = 1:length(subfolders)
        subfolderPath = fullfile(currentPath, subfolders{j});

        % Check for "_1b.mat" files in the current subfolder.
        fileList = dir(fullfile(subfolderPath, '*_1b.mat'));

        if isempty(fileList)
            fprintf('Folder "%s" does not contain a file ending with "_1b.mat".\n', subfolderPath);
        end
    end
end
end