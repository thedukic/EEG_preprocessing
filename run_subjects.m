function [errorLog, list_failed] = run_subjects(myPaths, errorLog)
% PROCESS_SUBJECT_BATCH Runs cleaning steps 1 and 2 for a batch of subjects,
% capturing errors into a centralized log without stopping execution.

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

    % Report
    report_final(myPaths, myPaths.subjects);
else
    warning('No participants selected: %s T%d', myPaths.group, myPaths.visit);
end

end