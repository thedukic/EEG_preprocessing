function log_error = run_test_method(myPaths)
% close all hidden; fclose all; clear all; clc;
% myPaths = preproc_folders;

% Initialise
log_error = struct('group', {}, 'visit', {}, 'subject', {}, 'index', [], 'step', {}, 'message', {});
path_plot = 'C:\DATA\MATLAB\myCodes\preprocessing\test\figures\test_ica';

% Select participants
i_group = 1;
i_visit = 1;
myPaths = preproc_participants(i_group, i_visit, myPaths);

% myPaths.subjects = {'ALS25976','ALS27315','ALS27392','ALS33970'};
% myPaths.subjects = {'ALS25976','ALS35439','ALS36062','ALS08518','ALS08585','ALS08637','ALS25511','ALS35439','ALS35973','ALS30088','ALS35862'};
% myPaths.subjects = {'ALS24900','ALS08707','ALS08549','ALS35225','ALS36431','ALS36062','ALS35996'};

% Sample size
num_subjects = length(myPaths.subjects);

for i_subj = 1:num_subjects
    % Header
    fprintf('\n');
    disp('==================================================================');
    disp([myPaths.task ' | ' ['T' num2str(myPaths.visit)] ' | ' myPaths.group ' | [' num2str(i_subj) '/' num2str(num_subjects) '] ' myPaths.subjects{i_subj} ' has started.']);
    disp('==================================================================');
    fprintf('\n');

    % Initialise
    close all hidden; fclose all; warning on;

    % Test
    try
        % run_test(myPaths, myPaths.subjects{i_subj}, path_plot);
        run_test_ica(myPaths, myPaths.subjects{i_subj}, path_plot);
    catch ME
        % Capture failure details
        newErr = struct(...
            'group', myPaths.group, ...
            'visit', ['T' num2str(myPaths.visit)], ...
            'subject', myPaths.subjects{i_subj}, ...
            'index', i_subj, ...
            'step', 'preproc_subject_1', ...
            'message', ME.message ...
            );
        log_error(end+1) = newErr; %#ok<AGROW>
        fprintf('Warning: Error encountered in cleaning step 1 for %s.\n', myPaths.subjects{i_subj});
    end
end


data_dir = 'C:\DATA\MATLAB\myCodes\preprocessing\test\figures\test_ica\classification';
[cohort_table, benchmark_stats] = benchmark_ica_cohort(data_dir);

end

% =========================================================================
% Helper functions
% =========================================================================
% function run_test(myPaths, id, path_plot)
%
% % Load preprocessing settings
% cfg = preproc_parameters;
%
% % -------------------------------------------------------------------------
% % Define paths and files
% % -------------------------------------------------------------------------
% subject = preproc_folders_subject(id, myPaths);
%
% % Print
% fprintf('==================================================================\n');
% fprintf('%s | %s | %s dataset | processing part 1 | pipeline v%s\n', subject.group, subject.id, subject.task, subject.codever);
% fprintf('==================================================================\n');
%
% % -------------------------------------------------------------------------
% % Basic data preparation
% % -------------------------------------------------------------------------
% % Find files
% subject.datablocks = list_datasets(subject.rawdata, subject.task);
%
% % Load those files
% EEG = load_biosemidata(subject, cfg);
% clearvars subject;
%
% % Fix events
% EEG = fix_events(EEG);
%
% % Add subject/channel info
% EEG = add_info(EEG);
%
% % Resample
% EEG = do_resampling(EEG, 256);
%
% % Keep event info
% EEG = extract_eventinfo(EEG, cfg.trg);
%
% % Filter highpass only
% EEG = do_filtering(EEG, 'highpass', cfg.flt);
%
% % Reference
% EEG = do_reref(EEG, 'aRobust');
%
% % -------------------------------------------------------------------------
% % Organise
% % -------------------------------------------------------------------------
% % Make EXT bipolar
% EEG = make_extbipolar(EEG);
%
% % Merge blocks
% EEG = merge_eeglabblocks(EEG);
%
% % Separate EXT
% [EEG, EMG, EXT] = separate_electrodetypes(EEG);
%
% % -------------------------------------------------------------------------
% % Test
% % -------------------------------------------------------------------------
% cfg_ecg = struct();
% cfg_ecg.win_heart    = [-200, 250];
% cfg_ecg.plot_visible = true;
% cfg_ecg.do_plot      = true;
% EXT.ALSUTRECHT.subject.figures = path_plot;
% detect_ecg_new(EXT, cfg_ecg);
%
% end

function run_test_ica(myPaths, id, path_plot)

% Time start
t_start = tic;

% Load preprocessing settings
cfg = preproc_parameters;
subject = preproc_folders_subject(id, myPaths);

% Modifications
subject.figures  = path_plot;
subject.data     = path_plot;
cfg.ica.num_ica  = 50;
cfg.flt.rs.hp(1) = 1;

% Print
fprintf('==================================================================\n');
fprintf('%s | %s | %s dataset | processing part 1 | pipeline v%s\n', subject.group, subject.id, subject.task, subject.codever);
fprintf('==================================================================\n');

% -------------------------------------------------------------------------
% Basic data preparation
% -------------------------------------------------------------------------
% % Find files
% subject.datablocks = list_datasets(subject.rawdata, subject.task);
%
% % Load those files
% EEG = load_biosemidata(subject, cfg);
% clearvars subject;
%
% % Fix events
% EEG = fix_events(EEG);
%
% % Add subject/channel info
% EEG = add_info(EEG);
%
% % Resample
% EEG = do_resampling(EEG, 256);
%
% % Keep event info
% EEG = extract_eventinfo(EEG, cfg);
%
% % Filter highpass only
% EEG = do_filtering(EEG, 'highpass', cfg);
%
% % Reference
% EEG = do_reref(EEG, 'aRobust');
%
% % Remove line noise
% EEG = reduce_linenoise1(EEG, cfg);
%
% % Make EXT bipolar
% EEG = make_extbipolar(EEG);
%
% % Merge blocks
% EEG = merge_eeglabblocks(EEG);
%
% % Separate EXT
% [EEG, EMG, EXT] = separate_electrodetypes(EEG);
%
% % Templates for ICA
% generate_ictemplateweights(EEG, EMG, EXT, cfg);
%
% % Remove noisy electrodes
% EEG = remove_noisyelectrodes(EEG, cfg);
%
% % Interpolate bad electrodes
% EEG = do_channelinterp(EEG, 'spherical');
%
% % Reference
% EEG = do_reref(EEG, 'aRegular');

% -------------------------------------------------------------------------
% Test 1
% -------------------------------------------------------------------------
% % ICA
% EEG = do_ica(EEG, cfg);
%
% % Save for next round!
% file_name = [EEG.ALSUTRECHT.subject.id '_' EEG.ALSUTRECHT.subject.visit '_' EEG.ALSUTRECHT.subject.task '_data_ica.mat'];
% save(fullfile(EEG.ALSUTRECHT.subject.data, file_name), "EEG", "EXT", "EMG");

% Load
file_name = [subject.id '_' subject.visit '_' subject.task '_data_ica.mat'];
load(fullfile(subject.data, 'data', file_name), "EEG", "EXT", "EMG");
EEG.ALSUTRECHT.subject.data = fullfile(subject.data, 'data');
EXT.ALSUTRECHT.subject.data = fullfile(subject.data, 'data');

% Detect artifact ICs
EEG2 = detect_badcomponents_new(EEG, EXT, EMG, cfg);

% % Remove artifact ICs
% EEG2 = remove_badcomponents(EEG2, cfg);
%
% % Report ICA
% report_ica(EEG2, cfg);

ica_class = EEG2.ALSUTRECHT.ica;
path_save = 'C:\DATA\MATLAB\myCodes\preprocessing\test\figures\test_ica\classification';
file_name = [EEG2.ALSUTRECHT.subject.id '_ica_class.mat'];
save(fullfile(path_save, file_name), "ica_class");

% -------------------------------------------------------------------------
% Test 2
% -------------------------------------------------------------------------
% % Merge
% EEG2 = merge_electrodetypes(EEG2, EMG, EXT);
%
% % MWF / DSS
% EEG2 = detect_blink_residuals(EEG2);
%
% if EEG2.ALSUTRECHT.leftovers.blink_info.has_residue
%     % EEG2 = do_mwf_blink(EEG2, cfg);
%     EEG2 = do_dss_blink(EEG2, cfg);
% else
%     fprintf('Skipping DSS as there was no evidence of blink leftovers.\n');
% end

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
elapsed_sec = toc(t_start);
if elapsed_sec < 60
    fprintf('\nElapsed time: %.2f seconds.\n', elapsed_sec);
else
    fprintf('\nElapsed time: %.2f seconds (%.2f minutes).\n', elapsed_sec, elapsed_sec / 60);
end

end