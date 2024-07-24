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

drivedata = 'E:';
myfolders.excpath = fullfile(drivedata,'2_OTHER_DATA\Excel\Utrecht\');
myfolders.gendata = [myfolders.excpath 'C9STATUS.xlsx'];
myfolders.peddata = [myfolders.excpath 'EEGPED.xlsx'];
myfolders.cogdata = [myfolders.excpath 'ECAS.xlsx'];
myfolders.nexdata = [myfolders.excpath 'NE.xlsx'];
myfolders.dmdata1 = [myfolders.excpath 'Table1_LME.xlsx'];
addpath('C:\DATA\MATLAB\myCodes\RS\common');

% Who to processed?
% 1. Selected or all participants
myfiles.todo = {};          % If empty, all are (re)done

% 2. C9 without the excluded participants
% % 1. Too noisy data
% % 2. Psychoactive medication
% % 3. Headtrauma
% excl_medication = {'ALS26360','ALS35895','ALS08665'};
% excl_headtrauma = {'ALS34168'}; % Does not have ECAS anyway ...
% cfg.excl        = [excl_medication,excl_headtrauma];
% myfiles.todo = select_participants(cfg,'C9',myfolders);

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
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
close all;