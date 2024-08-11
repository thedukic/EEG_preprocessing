
% =========================================================================
%
% EEG preprocessing main file, ALS Centre UMC Utrecht
% Check README.md for instructions
% SDukic, August 2024
%
% =========================================================================

close all; fclose all; clc; clearvars summaries;
delete(gcp("nocreate")); parpool("Processes");
myPaths = preproc_folders;

% Add specific info: group/task/visit
for i = 1  % :length(myPaths.group)
    for j = 1 % :length(myPaths.visit)
        myPathsTmp          = myPaths;
        myPathsTmp.task     = myPaths.task;
        myPathsTmp.group    = myPaths.group{i};
        myPathsTmp.visit    = myPaths.visit{j};
        myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata,myPathsTmp.group,myPathsTmp.visit);
        myPathsTmp.preproc  = fullfile(myPathsTmp.rootpreproc,myPathsTmp.task,myPathsTmp.group,myPathsTmp.visit);

        % Be careful when defining the to-do list
        % if processing more than one group
        subjects = list_subjects(myPathsTmp.rawdata,{});
        % myPathsTmp.excl = {'ALS08665','ALS34168'};
        % subjects = select_participants([],'C9',myPathsTmp);
        NSUB = length(subjects);

        if NSUB>0
            for k = 1:NSUB
                fprintf('\n');
                disp('==================================================================');
                disp([myPathsTmp.task ' | ' myPathsTmp.visit ' | ' myPathsTmp.group ' | [' num2str(k) '/' num2str(NSUB) '] ' subjects{k} ' has started.']);
                disp('==================================================================');
                fprintf('\n');

                % Clean 1&2
                output = preproc_cleaning1(myPathsTmp,subjects{k});
                output = preproc_cleaning2(myPathsTmp,subjects{k});

                % Record warnings for all participants in single table
                summaries(k,:) = struct2table(output,'AsArray',true);
            end

            % Report folder
            myPathsTmp.reports = fullfile(myPathsTmp.preproc,'reports');
            if exist(myPathsTmp.reports,'dir')~=7, mkdir(myPathsTmp.reports); end

            % Report table
            save(fullfile(myPathsTmp.reports,['Summary_' myPathsTmp.group '_' myPathsTmp.visit '_' myPathsTmp.task '_' myPathsTmp.proctime]),'summaries');
            writetable(summaries,fullfile(myPathsTmp.reports,['Summary_' myPathsTmp.group '_' myPathsTmp.visit '_' myPathsTmp.task  '_' myPathsTmp.proctime '.xlsx']),"WriteMode","overwrite");
            % clearvars summaries

            % Report visual
            report_final(myPathsTmp,1);
        end
    end
end

% =========================================================================