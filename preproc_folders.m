function myPaths = preproc_folders
% PREPROC_FOLDERS Configure paths, dependencies, cohort metadata, and parallel pools.
%
% Syntax:
%   myPaths = preproc_folders
%
% Description:
%   Initialises the execution environment for the automated EEG preprocessing
%   pipeline. Defines root input/output directories, sets cohort selection
%   parameters (tasks, study groups, longitudinal visits), resolves dependencies,
%   configures EEGLAB preferences, and manages parallel worker pools.
%
%   Key operations:
%     1. Sets root paths for pipeline source code, raw data, and preprocessed outputs.
%     2. Configures target cohort variables (task, groups, subgroups, visit range).
%     3. Recursively resolves and adds internal modules and third-party toolboxes
%        (EEGLAB, NoiseTools, Zapline-plus, GEDAI, restingIAF, BrewerMap) to MATLAB path.
%     4. Inspects path hierarchy for function collisions and duplicate definitions.
%     5. Validates storage drive existence and file-system accessibility.
%     6. Configures EEGLAB global memory and parallel computing preferences.
%     7. Purges orphaned cluster jobs/crash dumps and spins up a clean parallel pool.
%
% Inputs:
%   None (directory paths and cohort parameters are configured within this file).
%
% Outputs:
%   myPaths - Structure containing environment paths, cohort parameters, and metadata:
%               .codever     : Pipeline version string (e.g., '2')
%               .mycodes     : Root directory of pipeline source code
%               .rootrawdata : Input directory containing raw BioSemi recordings
%               .rootpreproc : Output directory for cleaned datasets and logs
%               .task        : Active task identifier (e.g., 'RS', 'MMN', 'SART', 'MT')
%               .group       : Cell array of target study groups (e.g., {'AFM'}, {'ALS'})
%               .subgroup    : Specific cohort subgroup identifier
%               .visit       : Numeric vector of session visits (e.g., 1:5)
%               .table1      : Path to exported clinical demographics table
%               .proctime    : Formatted timestamp string of pipeline launch
%
% ALS Centre, University Medical Centre Utrecht
% License: GNU General Public License v3.0

% =========================================================================
% !!!!!!!!!!!!      Define the parameters in this section      !!!!!!!!!!!!
% =========================================================================
% The scritop assumes that your data are strctured as, examples:
% E:\1_EEG_DATA\ALS\ALS12345\ALS12345_T1_MMN1.bdf
% E:\1_EEG_DATA\ALS\C42\C42_T1_EO1.bdf

% Root path of the pipeline code
myPaths.mycodes     = 'C:\DATA\MATLAB\myCodes\preprocessing';

% Root path containing raw EEG recordings
myPaths.rootrawdata = 'E:\1_EEG_DATA';

% Root path where preprocessed files will be saved
myPaths.rootpreproc = 'E:\3_PREPROCESSED_DATA';

% Target experimental paradigm (char): 'MMN', 'SART', 'RS', or 'MT'
myPaths.task        = 'MMN';

% Target cohort group(s) (cell array): {'ALS'}, {'CONTROL'}, {'AFM'}, {'PLS'}, or {'PMA'}
myPaths.group       = {'ALS'};

% Utrecht datasets: Optional subgroup filter (char): '', 'AFM_C9ORF72', 'AFM_ARPP21', 'MND_C9ORF72', or 'MND_SOD1'
myPaths.subgroup    = '';

% Longitudinal visit sessions to process (numeric vector): 1:5 or a single visit such as 1
myPaths.visit       = 1:5;

% Utrecht datasets: Master clinical metadata file exported from R
myPaths.table1      = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-08-14.txt';




% =========================================================================
% =========================================================================
% =========================================================================
% !!!!!!!!!!!!     The script below does not need changing     !!!!!!!!!!!!
% =========================================================================
% =========================================================================
% =========================================================================
warning on; warning('off', 'backtrace');

% Preprocessing code version
myPaths.codever = '2';

fprintf('==================================================================\n');
fprintf('Setting up the paths and loading the toolboxes\n');
fprintf('==================================================================\n');

fprintf('Pipeline version: %s\n\n', myPaths.codever);
fprintf('EEG data paths:\n');
fprintf('Raw: %s\n', myPaths.rootrawdata);
fprintf('Cleaned: %s\n', myPaths.rootpreproc);

% Track time
myPaths.proctime = strrep(strrep(char(datetime("now")), ':', '-'), ' ', '-');

% Navigate the main folder
cd(myPaths.mycodes);

% Clear exisiting paths
clean_toolbox_paths;

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
thisFolder      = pathsFolders{contains(pathsFolders, 'external', 'IgnoreCase', true)};
listFolders     = dir(thisFolder);
listFolders     = listFolders([listFolders.isdir]);
pathsFoldersTmp = fullfile(thisFolder, {listFolders(3:end).name});

if isempty(pathsFoldersTmp)
    thisFolder = 'C:\DATA\MATLAB\myCodes\external';
    fprintf('Your ''external'' folder is empty.\nUsing fallback instead: %s\n', thisFolder);

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
check_duplicate_all(myPaths.mycodes, {'external', 'unused'});
check_duplicate('preproc_main.m');
check_duplicate('run_subject_1.m');
check_duplicate('run_subject_2.m');
check_duplicate('eeglab.m');
check_duplicate('brewermap.m');
fprintf('\n');

% Initialise the toolboxes
eeglab; close all;

% Double-check drives
drive_1 = myPaths.rootrawdata(1:3);
drive_2 = myPaths.rootpreproc(1:3);
if ~(isfolder(drive_1) && isfolder(drive_2))
    error('Data paths are not correct. These local/online drives (%s or %s) do not exist.', drive_1, drive_2);
end

% Set EEGLAB options
pop_editoptions( ...
    'option_parallel', 1, ...
    'option_single', 0, ...
    'option_computeica', 0 ...
    );

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