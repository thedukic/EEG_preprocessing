function check_DUBfolders(folderPath)
% Lists all folders in the given path and checks file counts and naming.

if ~isfolder(folderPath)
    error('%s is not a valid folder.', folderPath);
end

folders = dir(folderPath);
folders = folders([folders.isdir] & ~strcmp({folders.name}, '.') & ~strcmp({folders.name}, '..')); % Filter for directories

for i = 1:length(folders)
    currentFolder = folders(i).name;
    currentFolderPath = fullfile(folderPath, currentFolder);
    files = dir(fullfile(currentFolderPath, '*'));
    files = files(~[files.isdir]); % Filter out subdirectories
    fileCount = length(files);

    fprintf('Folder: %s\n', currentFolder);
    if fileCount == 3
        fprintf('  Number of files: %d\n', fileCount);
    else
        warning('  Number of files: %d', fileCount);
    end

    for j = 1:fileCount
        currentFile = files(j).name;
        if ~isempty(currentFile)
            % Check the ID tags
            LID = length(currentFolder);
            subjID = currentFile(1:LID);
            if strcmpi(currentFolder, subjID)
                fprintf('  Good! Folder: %s -> File: %s \n', currentFile, currentFolder);
            else
                warning('  Mistake! Folder: %s -> File: %s \n', currentFile, currentFolder);
            end

            % Check the Tx tags
            Tx = currentFile(LID+2:LID+3);
            if ~strcmpi(Tx,'T1')
                warning('  Typo in the Tx tag? File: %s \n', currentFile);
            end

        else
            error('  Folder is empty???\n');
        end
    end
end

end