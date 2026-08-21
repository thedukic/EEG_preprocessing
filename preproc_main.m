% =========================================================================
%
% EEG preprocessing main file, ALS Centre UMC Utrecht
% Check README.md for instructions
% SDukic, July 2026
%
% TODO
% 1. Deal with cases where 1 file has 2 different tasks (MMN + SART)
% 2. https://github.com/bigdelys/eye-catch
% 3. Steamline the bad ic detection
% 4. Prevent outliers in EOG/ECG detection func

% =========================================================================
close all hidden; fclose all; clear all; clc;

% Initialise
myPaths = preproc_folders;
errorLog = struct('group', {}, 'visit', {}, 'subject', {}, 'index', [], 'step', {}, 'message', {});
listFailed = cell(length(myPaths.group), length(myPaths.visit));

% Loop
for i_group = 1:length(myPaths.group)
    for i_visit = 1:length(myPaths.visit)
        % Select participants
        myPathsTmp = preproc_participants(i_group, i_visit, myPaths);

        % Run
        [errorLog, listFailed{i_group, i_visit}] = run_subjects(myPathsTmp, errorLog);
    end
end

% Failure report
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
%
% generate_spatialdecay_profile(myPathsTmp, subjects);