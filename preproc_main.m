% =========================================================================
%
% EEG preprocessing main file, ALS Centre UMC Utrecht
% Check README.md for instructions
% SDukic, October 2024
%
% =========================================================================

close all; fclose all; clc; clear all;
myPaths = preproc_folders;

% Add specific info: group/task/visit
for i = 3     % :length(myPaths.group)
    for j = 1 % :length(myPaths.visit)
        myPathsTmp          = myPaths;
        myPathsTmp.task     = myPaths.task;
        myPathsTmp.group    = myPaths.group{i};
        myPathsTmp.visit    = myPaths.visit{j};
        myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata,myPathsTmp.group,myPathsTmp.visit);
        myPathsTmp.preproc  = fullfile(myPathsTmp.rootpreproc,myPathsTmp.task,myPathsTmp.group,myPathsTmp.visit);

        % Select participants
        [subjects, NSUB] = select_preproc_participants(myPathsTmp);

        if NSUB>0
            for k = 1:NSUB
                fprintf('\n');
                disp('==================================================================');
                disp([myPathsTmp.task ' | ' myPathsTmp.visit ' | ' myPathsTmp.group ' | [' num2str(k) '/' num2str(NSUB) '] ' subjects{k} ' has started.']);
                disp('==================================================================');
                fprintf('\n');

                % Cleaning step 1 & 2
                preproc_cleaning1(myPathsTmp,subjects{k});
                preproc_cleaning2(myPathsTmp,subjects{k});
            end

            % Report
            report_final(myPathsTmp,subjects,1);
            % estimate_ictemplates(myPathsTmp,subjects,1);
        end
    end
end
