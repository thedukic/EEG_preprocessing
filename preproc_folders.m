function [myfolders, myfiles] = preproc_folders

% Set paths
myfolders.mycodes     = 'C:\DATA\MATLAB\myCodes\SDukic\Preprocessing';
% myfolders.mytlbxs     = 'C:\DATA\MATLAB\Toolboxes';
myfolders.rootrawdata = 'E:\1_EEG_DATA';
myfolders.rootpreproc = 'E:\3_PREPROCESSED_DATA';

% Set group/task/visit
myfiles.task  = 'RS'; % MMN/SART/RS/MT
myfiles.group = {'AFM','ALS','PLS','PMA','CONTROL'};
myfiles.visit = {'T1'};   % T1/T2/...

% Add subfolders
files = dir(myfolders.mycodes);
subFolders = files([files.isdir]);
subFolderNames = {subFolders(3:end).name};
% subFolderPaths = strcat(myfolders.mycodes,subFolderNames);
subFolderPaths = [myfolders.mycodes, fullfile(myfolders.mycodes,subFolderNames)];
addpath(subFolderPaths{:});

% Add toolboxes from the external subfolder
subFolderExternal = subFolderPaths{contains(subFolderPaths,'external','IgnoreCase',true)};
files = dir(subFolderExternal);
subFolders = files([files.isdir]);
subFolderNames = {subFolders(3:end).name};
% subFolderPaths = strcat(subFolderExternal,subFolderNames);
subFolderPaths = fullfile(subFolderExternal,subFolderNames);
addpath(subFolderPaths{:});

if any(contains(subFolderPaths,'mwf'))
    addpath(genpath(subFolderPaths{contains(subFolderPaths,'mwf')}));
end

% addpath([myfolders.mytlbxs 'NoiseTools']);
% addpath([myfolders.mytlbxs 'eeglab2023.1']);
% addpath([myfolders.mytlbxs 'fieldtrip-20230913']);

% Initialise the toolboxes
% ft_defaults; ft_hastoolbox('cellfunction',1);
eeglab; close all;