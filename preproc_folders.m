function [myfolders, myfiles] = preproc_folders

% Set paths
myfolders.mycodes     = 'C:\DATA\MATLAB\myCodes\Preprocessing';                                                        % Pipeline
% myfolders.rootrawdata = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\1_RAW\EEG_DATA'; % Input
% myfolders.rootpreproc = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\2_PREPROCESSED'; % Output
myfolders.rootrawdata = 'E:\1_EEG_DATA'; % Input
myfolders.rootpreproc = 'E:\3_PREPROCESSED_DATA'; % Output

% Set group/task/visit
myfiles.task  = 'SART';          % MMN/SART/RS/EO/EC/MT
myfiles.group = {'ALS','CONTROL','AFM','PLS','PMA'};
myfiles.visit = {'T1','T2','T3','T4','T5'};

% Who to processed?
% myfiles.todo = list_missing('E:\1_EEG_DATA\AFM\T1','E:\3_PREPROCESSED_DATA\RS\AFM\T1');
% load('E:\4_POSTPROCESSED_DATA\AFM\EO\data\AFM_databag.mat','subjects');
% myfiles.todo = subjects;
myfiles.todo  = {'ALS26029'};          % If empty, all are (re)done

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
[ALLEEG, EEG, CURRENTSET] = eeglab;
close all;