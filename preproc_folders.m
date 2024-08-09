function myPaths = preproc_folders

% Track time
myPaths.proctime = strrep(strrep(char(datetime("now")),':','-'),' ','-');

% Set paths
myPaths.mycodes     = 'C:\DATA\MATLAB\myCodes\Preprocessing';                                                        % Pipeline
% myPaths.rootrawdata = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\1_RAW\EEG_DATA'; % Input
% myPaths.rootpreproc = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\2_PREPROCESSED'; % Output
myPaths.rootrawdata = 'D:\1_EEG_DATA';          % Input
myPaths.rootpreproc = 'D:\3_PREPROCESSED_DATA'; % Output

% Set task/group/visit
myPaths.task  = 'RS'; % MMN/SART/RS/EO/EC/MT
myPaths.group = {'ALS','CONTROL','AFM','PLS','PMA'};
myPaths.visit = {'T1','T2','T3','T4','T5'};

% drivedata = 'E:';
% myPaths.excpath = fullfile(drivedata,'2_OTHER_DATA\Excel\Utrecht\');
% myPaths.gendata = [myPaths.excpath 'C9STATUS.xlsx'];
% myPaths.peddata = [myPaths.excpath 'EEGPED.xlsx'];
% myPaths.cogdata = [myPaths.excpath 'ECAS.xlsx'];
% myPaths.nexdata = [myPaths.excpath 'NE.xlsx'];
% myPaths.dmdata1 = [myPaths.excpath 'Table1_LME.txt'];
% addpath('C:\DATA\MATLAB\myCodes\RS\common');

% Navigate the main folder
cd(myPaths.mycodes);

% Add subfolders
files             = dir(myPaths.mycodes);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = [myPaths.mycodes, fullfile(myPaths.mycodes,subFolderNames)];
addpath(subFolderPaths{:});

% Add toolboxes from the external subfolder
subFolderExternal = subFolderPaths{contains(subFolderPaths,'external','IgnoreCase',true)};
files             = dir(subFolderExternal);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = fullfile(subFolderExternal,subFolderNames);
addpath(subFolderPaths{:});

% Add MWF subfolders as well
if any(contains(subFolderPaths,'mwf'))
    addpath(genpath(subFolderPaths{contains(subFolderPaths,'mwf')}));
end

% Initialise the toolboxes
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; close all;

end