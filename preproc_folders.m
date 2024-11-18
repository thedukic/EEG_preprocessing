function myPaths = preproc_folders

% Track time
myPaths.proctime = strrep(strrep(char(datetime("now")),':','-'),' ','-');

% Set paths
myPaths.mycodes     = 'C:\DATA\MATLAB\myCodes\Preprocessing';                                                        % Pipeline
% myPaths.rootrawdata = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\1_RAW\EEG_DATA'; % Input
% myPaths.rootpreproc = 'L:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\2_PREPROCESSED'; % Output
myPaths.rootrawdata = 'E:\1_EEG_DATA';          % Input
myPaths.rootpreproc = 'E:\3_PREPROCESSED_DATA'; % Output

% Set task/group/visit
myPaths.task  = 'SART'; % MMN/SART/RS/EO/EC/MT
myPaths.group = {'ALS','CONTROL','AFM','PLS','PMA'};
myPaths.visit = {'T1','T2','T3','T4','T5'};

drivedata = 'E:';
myPaths.excpath = fullfile(drivedata,'2_OTHER_DATA\Excel\Utrecht\');
myPaths.gendata = [myPaths.excpath 'C9status.xlsx'];
myPaths.peddata = [myPaths.excpath 'Pedigrees.xlsx'];
myPaths.cogdata = [myPaths.excpath 'ECAS.txt'];
myPaths.nexdata = [myPaths.excpath 'NE.txt'];
myPaths.dmdata1 = [myPaths.excpath 'Table1.txt'];

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

% % Add MWF subfolders as well
% if any(contains(subFolderPaths,'mwf'))
%     addpath(genpath(subFolderPaths{contains(subFolderPaths,'mwf')}));
% end

% Initialise the toolboxes
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; close all;

% Double-check drives
drive1 = myPaths.rootrawdata(1:3);
drive2 = myPaths.rootpreproc(1:3);
if ~(isfolder(drive1) && isfolder(drive2))
    error('Data paths are not correct. These local/online drives (%s or %s) do not exist.',drive1,drive2);
end

% Set EEGLAB options
if strcmpi(myPaths.task,'MT')
    % Motor task data are large, and the current implementation of parallel processing in EEGLAB would require a large amount of RAM
    flagParallel = 0;
else
    flagParallel = 1;
end
pop_editoptions('option_parallel',flagParallel,'option_single',0,'option_computeica',0);

% Start parallel processes
pool = gcp('nocreate');
delete(pool); parpool("Processes");

% if ~isempty(pool)
%     if ~isempty(pool.Cluster) % pool.Cluster.HasSharedFilesystem && pool.SpmdEnabled
%         % disp('Running with processes.');
%     else
%         % disp('Running with threads.');
%         delete(pool); parpool("Processes");
%     end
% else
%     % disp('No active parpool.');
%     parpool("Processes");
% end

end