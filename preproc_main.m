% =========================================================================
%
% EEG preprocessing main file, ALS Centre UMC Utrecht
% Check README.md for instructions
% SDukic, August 2024
%
% =========================================================================

close all; fclose all; clc; clear all;
myPaths = preproc_folders;

delete(gcp("nocreate")); parpool("Processes");
pop_editoptions('option_parallel',1,'option_single',0);

% Add specific info: group/task/visit
for i = 3  % :length(myPaths.group)
    for j = 1 % :length(myPaths.visit)
        myPathsTmp          = myPaths;
        myPathsTmp.task     = myPaths.task;
        myPathsTmp.group    = myPaths.group{i};
        myPathsTmp.visit    = myPaths.visit{j};
        myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata,myPathsTmp.group,myPathsTmp.visit);
        myPathsTmp.preproc  = fullfile(myPathsTmp.rootpreproc,myPathsTmp.task,myPathsTmp.group,myPathsTmp.visit);

        % Be careful when defining the to-do list if processing more than one group!
        % myPathsTmp.excl = {}; subjects = select_participants([],'C9',myPathsTmp);
        subjects = list_subjects(myPathsTmp.rawdata,{});
        % subjects = list_subjects('E:\3_PREPROCESSED_DATA\RS\ALS\T1',{});
        % subjects = {'ALS26603','ALS37919'};
        NSUB = length(subjects);

        fprintf('Processing %d %s participants....\n',NSUB,myPathsTmp.group);
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
                if ~isempty(output)
                    summaries(k,:) = struct2table(output,'AsArray',true);
                else
                    summaries(end+1,:) = cell2table([subjects{k}, repmat({NaN},1,length(summaries(1,:).Properties.VariableNames)-1)], 'VariableNames', summaries(1,:).Properties.VariableNames);
                end
            end

            % Report folder
            myPathsTmp.reports = fullfile(myPathsTmp.preproc,'reports');
            if exist(myPathsTmp.reports,'dir')~=7, mkdir(myPathsTmp.reports); end

            % Report table
            save(fullfile(myPathsTmp.reports,['Summary_' myPathsTmp.group '_' myPathsTmp.visit '_' myPathsTmp.task '_' myPathsTmp.proctime]),'summaries');
            writetable(summaries,fullfile(myPathsTmp.reports,['Summary_' myPathsTmp.group '_' myPathsTmp.visit '_' myPathsTmp.task  '_' myPathsTmp.proctime '.xlsx']),"WriteMode","overwrite");
            clearvars summaries;

            % Report visual
            report_final(myPathsTmp,1);
        end
    end
end

% =========================================================================