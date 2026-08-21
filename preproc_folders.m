function myPaths = preproc_folders
% =========================================================================
%
% Script for setting up the paths and labels of data for preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================

% Preprocessing code version
myPaths.codever = '2';

% Define
myPaths.mycodes     = 'C:\DATA\MATLAB\myCodes\preprocessing';     % Pipeline
myPaths.rootrawdata = 'E:\1_EEG_DATA';                            % Input
myPaths.rootpreproc = 'E:\3_PREPROCESSED_DATA';                   % Output
% myPaths.rootrawdata = 'C:\DATA\MATLAB\EEG\1_EEG_DATA';          % Input
% myPaths.rootpreproc = 'C:\DATA\MATLAB\EEG\3_PREPROCESSED_DATA'; % Output

% Task (char): MMN / SART / RS / MT
myPaths.task  = 'RS';
% Group (cell): ALS / CONTROL / AFM / PLS / PMA
myPaths.group = {'AFM'};
% Subgroup (char): AFM_C9ORF72 / AFM_ARPP21 / MND_C9ORF72 / MND_SOD1
myPaths.subgroup = 'AFM_C9ORF72';
% Visit (num): 1-5
myPaths.visit = 1:5;

% Path to the table 1 (export from R)
myPaths.table1 = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-08-14.txt';

% =========================================================================
% The script below does not need changing
% =========================================================================
warning on; warning('off', 'backtrace');

fprintf('==================================================================\n');
fprintf('Setting up the paths and loading the toolboxes\n');
fprintf('==================================================================\n');
fprintf('Pipeline version: %s\n\n', myPaths.codever);

% Track time
myPaths.proctime = strrep(strrep(char(datetime("now")), ':', '-'), ' ', '-');

% Navigate the main folder
cd(myPaths.mycodes);

fprintf('EEG data paths:\n');
fprintf('Raw: %s\n', myPaths.rootrawdata);
fprintf('Cleaned: %s\n', myPaths.rootpreproc);

% Add folders
listFolders  = dir(myPaths.mycodes);
listFolders  = listFolders([listFolders.isdir]);
listFolders  = {listFolders(3:end).name};
pathsFolders = [myPaths.mycodes, fullfile(myPaths.mycodes, listFolders)];
pathsFolders = pathsFolders(~contains(pathsFolders, {'git', 'unused'}));

addpath(pathsFolders{:});
fprintf('Adding folders:\n');
fprintf('%s\n', pathsFolders{:});

% Add toolboxes from the "external" folder
thisFolder      = pathsFolders{contains(pathsFolders,'external','IgnoreCase',true)};
listFolders     = dir(thisFolder);
listFolders     = listFolders([listFolders.isdir]);
pathsFoldersTmp = fullfile(thisFolder, {listFolders(3:end).name});

if isempty(pathsFoldersTmp)
    thisFolder = 'C:\DATA\MATLAB\myCodes\external';
    fprintf('Your ''external'' folder is empty.\nUsing instead: %s\n', thisFolder);

    pathsFoldersTmp    = {};
    pathsFoldersTmp{1} = fullfile(thisFolder, 'eeglab2025.1.0');
    pathsFoldersTmp{2} = fullfile(thisFolder, 'noisetools_29-Apr-2023');
    pathsFoldersTmp{3} = fullfile(thisFolder, 'zaplineplus_14-Apr-2023');
    pathsFoldersTmp{4} = fullfile(thisFolder, 'gedai_05082026');
    pathsFoldersTmp{5} = fullfile(thisFolder, 'restingiaf_20-Jan-2025');
    pathsFoldersTmp{6} = fullfile(thisFolder, 'brewermap-3.2.8');
end

addpath(pathsFoldersTmp{:});
fprintf('Adding external toolboxes:\n');
fprintf('%s\n', pathsFoldersTmp{:});

% Add subfolers from the "files" folder
thisFolder      = pathsFolders{contains(pathsFolders,'files','IgnoreCase',true)};
listFolders     = dir(thisFolder);
listFolders     = listFolders([listFolders.isdir]);
pathsFoldersTmp = fullfile(thisFolder,{listFolders(3:end).name});

addpath(pathsFoldersTmp{:});
fprintf('Adding ''file'' subfolders:\n');
fprintf('%s\n', pathsFoldersTmp{:});

% Check for duplicates to prevent overloading
% restoredefaultpath % Maybe better not to use it altough it does the job
% Check if there are functions with the same name
check_duplicates(myPaths.mycodes, {'external','unused'});
check_duplicatefunc('preproc_main.m');
check_duplicatefunc('preproc_cleaning1.m');
check_duplicatefunc('preproc_cleaning2.m');
check_duplicatefunc('eeglab.m');
check_duplicatefunc('brewermap.m');
fprintf('\n');

% Initialise the toolboxes
eeglab; close all;

% Double-check drives
drive1 = myPaths.rootrawdata(1:3);
drive2 = myPaths.rootpreproc(1:3);
if ~(isfolder(drive1) && isfolder(drive2))
    error('Data paths are not correct. These local/online drives (%s or %s) do not exist.',drive1,drive2);
end

% Set EEGLAB options
pop_editoptions( ...
    'option_parallel', 1, ...
    'option_single', 0, ...
    'option_computeica',0);

% -------------------------------------------------------------------------
% Reset parallel architecture and clear legacy crash dumps
% -------------------------------------------------------------------------
% 1. Shut down any active or hanging pool first to release file locks
existingPool = gcp('nocreate');
if ~isempty(existingPool)
    delete(existingPool);
end

% 2. Access the cluster profile to clean up the workspace disk cache
myCluster = parcluster('Processes');
crashedJobs = myCluster.Jobs;

if ~isempty(crashedJobs)
    try
        delete(crashedJobs);
        fprintf('Successfully cleared %d legacy crash dump directories.\n', length(crashedJobs));
    catch
        % Guard against rare OS file-system locking delays
        warning('Some crash logs are currently locked by the OS and will be cleared next run.');
    end
end

% 3. Spin up a fresh, clean parallel pool
parpool("Processes");

end