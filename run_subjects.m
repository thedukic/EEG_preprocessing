function [errorLog, list_failed] = run_subjects(myPaths, errorLog)
% RUN_SUBJECTS Batch runner for subject-level preprocessing and postprocessing.
%
% Syntax:
%   [errorLog, list_failed] = run_subjects(myPaths, errorLog)
%
% Description:
%   Executes Part 1 (continuous cleaning) and Part 2 (postprocessing and QA)
%   sequentially for each participant in a selected cohort batch.
%
%   Key operations:
%     1. Iterates across all participant IDs defined in myPaths.subjects.
%     2. Wraps run_subject_1 and run_subject_2 in independent try-catch blocks
%        to prevent unhandled runtime exceptions from halting batch execution.
%     3. Logs failure metadata (group, visit, subject ID, step, error message)
%        into a centralised structured error log.
%     4. Automatically skips Part 2 if Part 1 fails for a given participant.
%     5. Evaluates output file completion via check_runs and compiles batch-level
%        summary reports via report_final.
%
% Inputs:
%   myPaths  - Structure containing batch and cohort metadata:
%                .subjects    : Cell array of participant ID strings
%                .group       : Study cohort or group identifier
%                .visit       : Target visit number (scalar, e.g., 1, 2)
%                .task        : Task identifier (e.g., 'RS', 'EO', 'EC', 'MT')
%                .rootrawdata : Root directory containing raw recordings
%                .rootpreproc : Root directory for preprocessed outputs
%                .mycodes     : Root directory of pipeline source code
%   errorLog - Struct array tracking processing errors across batch runs.
%
% Outputs:
%   errorLog    - Updated struct array with newly appended failure records.
%   list_failed - List of participants/sessions identified with incomplete runs.
%
% ALS Centre, University Medical Centre Utrecht
% License: GNU General Public License v3.0

% Sample size
num_subjects = length(myPaths.subjects);
list_failed = [];

if num_subjects > 0
    for i_subj = 1:num_subjects
        % Header
        fprintf('\n');
        disp('==================================================================');
        disp([myPaths.task ' | ' ['T' num2str(myPaths.visit)] ' | ' myPaths.group ' | [' num2str(i_subj) '/' num2str(num_subjects) '] ' myPaths.subjects{i_subj} ' has started.']);
        disp('==================================================================');
        fprintf('\n');

        % Initialise
        step1_passed = true;
        close all hidden; fclose all; warning on;

        % --- PART 1 ---
        try
            run_subject_1(myPaths, myPaths.subjects{i_subj});
        catch ME
            step1_passed = false;

            % Capture failure details
            newErr = struct(...
                'group', myPaths.group, ...
                'visit', ['T' num2str(myPaths.visit)], ...
                'subject', myPaths.subjects{i_subj}, ...
                'index', i_subj, ...
                'step', 'preproc_subject_1', ...
                'message', ME.message ...
                );
            errorLog(end+1) = newErr; %#ok<AGROW>
            fprintf('Warning: Error encountered in cleaning step 1 for %s.\n', myPaths.subjects{i_subj});
        end

        % --- PART 2 ---
        if step1_passed
            try
                run_subject_2(myPaths, myPaths.subjects{i_subj});
            catch ME
                % Capture failure details
                newErr = struct(...
                    'group', myPaths.group, ...
                    'visit', ['T' num2str(myPaths.visit)], ...
                    'subject', myPaths.subjects{i_subj}, ...
                    'index', i_subj, ...
                    'step', 'preproc_subject_2', ...
                    'message', ME.message ...
                    );
                errorLog(end+1) = newErr; %#ok<AGROW>
                fprintf('Warning: Error encountered in cleaning step 2 for %s.\n', myPaths.subjects{i_subj});
            end
        else
            warning on;
            fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            fprintf('Warning: Skipping cleaning step 2 because step 1 failed for %s.\n', myPaths.subjects{i_subj});
            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        end
    end

    % Check if any runs failed
    list_failed = check_runs(myPaths);

    % Report (you want a whole cohort)
    if num_subjects > 1
        report_final(myPaths, myPaths.subjects);
    end
else
    warning('No participants selected: %s T%d', myPaths.group, myPaths.visit);
end

end