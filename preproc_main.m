% =========================================================================
% PREPROC_MAIN Batch execution script for the automated EEG preprocessing pipeline.
%
% Syntax:
%   preproc_main
%
% Description:
%   Runs batch execution of the two-stage EEG preprocessing pipeline
%   across cohorts, study groups, longitudinal visits, and tasks.
%
%   Key operations:
%     1. Environment and path initialisation via preproc_folders.
%     2. Participant selection via preproc_participants. 
%     3. Batch execution via run_subjects.
%     4. Structured error tracking and generation of failure logs.
%
% TODO:
%   1. Evaluate Eye-Catch integration for automated ocular IC detection.
%   2. Offsets not correct for dataset collected with fs > 256 Hz.
%
% ALS Centre, University Medical Centre Utrecht
% Author: S. Dukic, August 2026
% License: GNU General Public License v3.0
% =========================================================================

% =========================================================================
% 0. Pipeline setup via these files
% Read: EEG Preprocessing Pipeline Manual.md
% =========================================================================
% - preproc_folders.m
% - preproc_participants.m

% =========================================================================
% 1. Batch process: group & visit
% =========================================================================
% Housekeeping 
close all hidden; fclose all; clear all; clc;

% Initialise
myPaths    = preproc_folders;
errorLog   = struct('group', {}, 'visit', {}, 'subject', {}, 'index', [], 'step', {}, 'message', {});
listFailed = cell(length(myPaths.group), length(myPaths.visit));

% =========================================================================
% 2. Batch process: group & visit
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
% 2. Batch process: task (for pipeline testing)
% =========================================================================
% % Select participants
% list_tasks = {'MMN', 'SART', 'RS', 'MT'};
% i_group = 1;
% i_visit = 1;
% 
% % Run
% for i_task = 1:length(list_tasks)
%     % Select task
%     myPaths.task = list_tasks{i_task};
%     myPathsTmp = preproc_participants(i_group, i_visit, myPaths);
% 
%     % Run
%     [errorLog, listFailed{i_task, i_visit}] = run_subjects(myPathsTmp, errorLog);
% end

% =========================================================================
% 3. Failure report
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