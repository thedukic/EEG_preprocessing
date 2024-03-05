%
% EEG preprocessing main file
%
% 1. CUDICA must be set up, see instructions in its folder
% 2. Set/check paths (and settings):
%    a. preproc_folders;
%    b. preproc_parameters;

%% ========================================================================
% Run the code below
[myfolders, myfiles] = preproc_folders;

for i = 1 % :length(myfiles.group)
    % Add specific info: group/task/visit
    myfolders.group   = myfiles.group{i};
    myfolders.task    = myfiles.task;
    myfolders.visit   = myfiles.visit{1};
    myfolders.rawdata = fullfile(myfolders.rootrawdata,myfolders.group,myfolders.visit);
    myfolders.preproc = fullfile(myfolders.rootpreproc,myfolders.task,myfolders.group,myfolders.visit);

    subjects = list_subjects(myfolders.rawdata,[]);
    for j = 6:length(subjects)
        output = preproc_cleaning(myfolders,subjects{j});

        % Record warnings for all participants in single table
        if isempty(output)
            summaries(j,1) = subjects(j);
            summaries(j,2:end) = {NaN};
        else
            summaries(j,:) = struct2table(output,'AsArray',true);
        end
    end
    save(fullfile(myfolders.preproc,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task]),'summaries');
    writetable(summaries,fullfile(myfolders.preproc,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task '.xlsx']));
end

%% ========================================================================