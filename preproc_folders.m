function myPaths = preproc_folders
% =========================================================================
%
% Script for setting up the paths and labels of data for preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================

% Preprocessing code version
myPaths.rnum = '1';

% Define
myPaths.mycodes     = 'C:\DATA\MATLAB\myCodes\Preprocessing';   % Pipeline
myPaths.rootrawdata = 'D:\1_EEG_DATA';          % Input
myPaths.rootpreproc = 'D:\3_PREPROCESSED_DATA'; % Output

% Task: MMN/SART/RS/MT
myPaths.task  = 'RS';
% Group: ALS / CONTROL / AFM / PLS/ PMA
myPaths.group = {'ALS'};
% Visit: T1/ T2 / T3 / T4 / T5
myPaths.visit = {'T1'};

% =========================================================================
% The script below does not need changing
% =========================================================================
% set(0, 'DefaultFigureVisible', 'off'); % not tested yet
warning on; warning('off','backtrace');

fprintf('\n');
disp('==================================================================');
fprintf('Setting up the paths and loading the toolboxes...\n');
disp('==================================================================');

% Track time
myPaths.proctime = strrep(strrep(char(datetime("now")),':','-'),' ','-');

% Navigate the main folder
cd(myPaths.mycodes);

fprintf('EEG data paths:\n');
fprintf('Raw: %s\n', myPaths.rootrawdata);
fprintf('Cleaned: %s\n', myPaths.rootpreproc);

% Add subfolders
files             = dir(myPaths.mycodes);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = [myPaths.mycodes, fullfile(myPaths.mycodes,subFolderNames)];
subFolderPaths    = subFolderPaths(~contains(subFolderPaths,{'git','unused'}));

addpath(subFolderPaths{:});
fprintf('Adding folders:\n');
fprintf('%s\n', subFolderPaths{:});

% Add toolboxes from the external subfolder
subFolderExternal = subFolderPaths{contains(subFolderPaths,'external','IgnoreCase',true)};
files             = dir(subFolderExternal);
subFolders        = files([files.isdir]);
subFolderNames    = {subFolders(3:end).name};
subFolderPaths    = fullfile(subFolderExternal,subFolderNames);

addpath(subFolderPaths{:});
fprintf('Adding external toolboxes:\n');
fprintf('%s\n', subFolderPaths{:}); fprintf('\n');

% Check for duplicates to prevent overloading
% restoredefaultpath % Maybe better not to use it altough it does the job
check_duplicatefunc('preproc_main.m');
check_duplicatefunc('preproc_cleaning1.m');
check_duplicatefunc('preproc_cleaning2.m');
check_duplicatefunc('eeglab.m');
check_duplicatefunc('brewermap.m');

% Initialise the toolboxes
eeglab; close all;

% Double-check drives
drive1 = myPaths.rootrawdata(1:3);
drive2 = myPaths.rootpreproc(1:3);
if ~(isfolder(drive1) && isfolder(drive2))
    error('Data paths are not correct. These local/online drives (%s or %s) do not exist.',drive1,drive2);
end

% Set EEGLAB options
if strcmpi(myPaths.task,'MT')
    % Motor task data are large
    % The current implementation of parallel processing in EEGLAB would require large RAM
    flagParallel = 0;
else
    flagParallel = 1;
end
pop_editoptions('option_parallel',flagParallel,'option_single',0,'option_computeica',0);

% Kill and start again the parallel processes
delete(gcp('nocreate')); parpool("Processes");

% The code is supposed to be smart about killing/starting the parallel processes
% But cannot make it work
% pool = gcp('nocreate');
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