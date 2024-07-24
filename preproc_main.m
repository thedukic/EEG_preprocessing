
% =========================================================================
%
% EEG preprocessing main file
% Check README.md for instructions
% SDukic, April 2024
%
% =========================================================================

close all; fclose all; clc; clearvars summaries;
[myfolders, myfiles] = preproc_folders;

proctime = strrep(strrep(char(datetime("now")),':','-'),' ','-');

% Add specific info: group/task/visit
for i = 3  % :length(myfiles.group)
    for j = 1 % :length(myfiles.visit)
        myfolders.proctime = proctime;
        myfolders.group    = myfiles.group{i};
        myfolders.task     = myfiles.task;
        myfolders.visit    = myfiles.visit{j};
        myfolders.rawdata  = fullfile(myfolders.rootrawdata,myfolders.group,myfolders.visit);
        myfolders.preproc  = fullfile(myfolders.rootpreproc,myfolders.task,myfolders.group,myfolders.visit);

        % To-do list makes sense if you set one cohort only!
        subjects = list_subjects(myfolders.rawdata,myfiles.todo);
        NSUB = length(subjects);

        if NSUB>0
            for k = 1:NSUB
                fprintf('\n');
                disp('==================================================================');
                disp([myfolders.task ' | ' myfolders.visit ' | ' myfolders.group ' | [' num2str(k) '/' num2str(NSUB) '] ' subjects{k} ' has started.']);
                disp('==================================================================');
                fprintf('\n');
                output = preproc_cleaning(myfolders,subjects{k});

                % Record warnings for all participants in single table
                summaries(k,:) = struct2table(output,'AsArray',true);
            end

            % Report
            myfolders.reports = fullfile(myfolders.preproc,'reports');
            if exist(myfolders.reports,'dir')~=7, mkdir(myfolders.reports); end

            % Visual report
            report_final(myfolders,1);

            % Add that there is a group-level overview of channels being rejected
            % Save the summary report
            save(fullfile(myfolders.reports,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task '_' myfolders.proctime]),'summaries');
            writetable(summaries,fullfile(myfolders.reports,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime '.xlsx']));
            clearvars summaries
        end
    end
end

% =========================================================================