function [myfolders, myfiles] = preproc_folders

% Set paths
myfolders.mycodes     = 'C:\DATA\MATLAB\myCodes\SDukic\Preprocessing';  % Pipeline
myfolders.rootrawdata = 'E:\1_EEG_DATA';                                % Input
myfolders.rootpreproc = 'E:\3_PREPROCESSED_DATA';                       % Output

% Set group/task/visit
myfiles.task  = 'RS';          % MMN/SART/RS/MT
myfiles.group = {'ALS','CONTROL','AFM','PLS','PMA'};
myfiles.visit = {'T1'};        % T1/T2/...
myfiles.todo  = {'ALS34543'};  % If empty, all are done

% Navigate the main folder
cd(myfolders.mycodes);

% Add subfolders
files             = dir(myfolders.mycodes);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = [myfolders.mycodes, fullfile(myfolders.mycodes,subFolderNames)];
addpath(subFolderPaths{:});

% Add toolboxes from the external subfolder
subFolderExternal = subFolderPaths{contains(subFolderPaths,'external','IgnoreCase',true)};
files             = dir(subFolderExternal);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = fullfile(subFolderExternal,subFolderNames);
addpath(subFolderPaths{:});

% Add MWF subfolders 
if any(contains(subFolderPaths,'mwf'))
    addpath(genpath(subFolderPaths{contains(subFolderPaths,'mwf')}));
end

% Initialise the toolboxes
eeglab; close all;