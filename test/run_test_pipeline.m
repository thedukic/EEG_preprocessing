function log_error = run_test_pipeline(myPaths)
% close all hidden; fclose all; clear all; clc;
% myPaths = preproc_folders;

% Select participants
list_tasks = {'MMN', 'SART', 'RS', 'MT'};
i_group = 1;
i_visit = 1;

% Initialise
log_error = struct('group', {}, 'visit', {}, 'subject', {}, 'index', [], 'step', {}, 'message', {});
list_fail = cell(length(list_tasks));

% Batch process: task
for i_task = 1:length(list_tasks)
    % Select task / participants
    myPaths.task = list_tasks{i_task};
    myPathsTmp = preproc_participants(i_group, i_visit, myPaths);

    % Run
    [log_error, list_fail{i_task}] = run_subjects(myPathsTmp, log_error);
end

end
