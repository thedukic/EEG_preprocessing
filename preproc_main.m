% =========================================================================
% PREPROC_MAIN Batch execution script for the automated EEG preprocessing pipeline.
%
% Syntax:
%   preproc_main
%
% Description:
%   Orchestrates batch execution of the two-stage EEG preprocessing pipeline
%   across cohorts, study groups, longitudinal visits, and tasks.
%
%   Key operations:
%     1. Environment initialisation and path resolution via preproc_folders.
%     2. Cohort iteration across defined study groups and visit sessions.
%     3. Participant selection and batch execution via run_subjects.
%     4. Structured error tracking and generation of overnight failure logs (CSV).
%
% Requirements:
%   - Paths must be configured in preproc_folders.m.
%   - Pipeline parameters must be defined in preproc_parameters.m.
%
% Outputs:
%   - Preprocessed and epoched EEG datasets (.set / .mat).
%   - Diagnostic figures and QA reports per subject/block.
%   - 'overnight_pipeline_errors.csv' written to rootpreproc upon completion.
%
% TODO:
%   1. Evaluate Eye-Catch integration for automated ocular IC detection.
%   2. Offsets not correct for dataset collected with fs > 256 Hz
%
% ALS Centre, University Medical Centre Utrecht
% Author: S. Dukic, August 2026
% License: GNU General Public License v3.0
% =========================================================================

close all hidden; fclose all; clear all; clc;

% Initialise
myPaths = preproc_folders;
errorLog = struct('group', {}, 'visit', {}, 'subject', {}, 'index', [], 'step', {}, 'message', {});
listFailed = cell(length(myPaths.group), length(myPaths.visit));

% =========================================================================
% Batch process: group & visit
% =========================================================================
for i_group = 1:length(myPaths.group)
    for i_visit = 1:length(myPaths.visit)
        % Select participants
        myPathsTmp = preproc_participants(i_group, i_visit, myPaths);

        % Run
        [errorLog, listFailed{i_group, i_visit}] = run_subjects(myPathsTmp, errorLog);
    end
end

% =========================================================================
% Batch process: task (for pipiline testing)
% =========================================================================
% Select participants
list_tasks = {'MMN', 'SART', 'RS', 'MT'};
i_group = 1;
i_visit = 1;

% Run
for i_task = 1:length(list_tasks)
    % Select task
    myPaths.task = list_tasks{i_task};
    myPathsTmp = preproc_participants(i_group, i_visit, myPaths);

    % Run
    [errorLog, listFailed{i_task, i_visit}] = run_subjects(myPathsTmp, errorLog);
end

% =========================================================================
% Failure report
% =========================================================================
fprintf('\n\n');
disp('==================================================================');
disp('PROCESSING FINISHED: FAILURE SUMMARY');
disp('==================================================================');
if isempty(errorLog)
    disp('All participants and processing steps completed successfully!');
else
    errorTable = struct2table(errorLog);
    disp(errorTable);

    writetable(errorTable, fullfile(myPaths.rootpreproc, 'overnight_pipeline_errors.csv'));
    fprintf('The failure log has been saved to overnight_pipeline_errors.csv\n');
end









% =========================================================================
% myPathsTmp.preproc = cell(1,2);
% myPathsTmp.preproc{1} = 'E:\3_PREPROCESSED_DATA\RS\CONTROL\T1';
% myPathsTmp.preproc{2} = 'E:\3_PREPROCESSED_DATA\RS\ALS\T1';
% report_final(myPaths, 'E:\3_PREPROCESSED_DATA\RS');
% generate_spatialdecay_profile(myPathsTmp, subjects);
% =========================================================================