
% =========================================================================
%
% EEG preprocessing main file
% Check README.md for instructions
% SDukic, April 2024
%
% =========================================================================

close all; fclose all; clc;
[myfolders, myfiles] = preproc_folders;

proctime = strrep(strrep(char(datetime("now")),':','-'),' ','-');

for i = 3 % :length(myfiles.group)
    % Add specific info: group/task/visit
    myfolders.group   = myfiles.group{i};
    myfolders.task    = myfiles.task;
    myfolders.visit   = myfiles.visit{1};
    myfolders.rawdata = fullfile(myfolders.rootrawdata,myfolders.group,myfolders.visit);
    myfolders.preproc = fullfile(myfolders.rootpreproc,myfolders.task,myfolders.group,myfolders.visit);

    % To-do list makes sense if you set one cohort only!
    subjects = list_subjects(myfolders.rawdata,myfiles.todo);
    NSUB = length(subjects);

    for j = 1:NSUB
        fprintf('\n'); 
        disp('==================================================================');
        disp([myfolders.task ' | ' myfolders.group ' | [' num2str(j) '/' num2str(NSUB) '] ' subjects{j} ' has started.']);
        disp('=================================================================='); 
        fprintf('\n'); 
        output = preproc_cleaning(myfolders,subjects{j});

        % Record warnings for all participants in single table
        summaries(j,:) = struct2table(output,'AsArray',true);
    end

    % Add that there is a group-level overview of channels being rejected
    % Save the summary report
    save(fullfile(myfolders.preproc,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task '_' proctime]),'summaries');
    writetable(summaries,fullfile(myfolders.preproc,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' proctime '.xlsx']));
    clearvars summaries
end

% =========================================================================