% =========================================================================
%
% EEG preprocessing main file, ALS Centre UMC Utrecht
% Check README.md for instructions
% SDukic, March 2025
%
% TODO
% 1. Turn off figure visibility while plotting (speed up the code)
% 2. ZipLine that exlcudes other peaks basides the line noise
% 3. When checking leftovers in step1, interpolate the outlier channels
% 4. Deal with files that have diff tasks in them, like MMN+SART
% 5.
% =========================================================================

close all; fclose all; clear all; clc;
myPaths = preproc_folders;

% Run
for i = 1:length(myPaths.group)
    for j = 1:length(myPaths.visit)
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

                % Cleaning steps
                preproc_cleaning1(myPathsTmp,subjects{k});
                preproc_cleaning2(myPathsTmp,subjects{k});
            end

            % Report
            report_final(myPathsTmp,subjects);
        end
    end
end
