function run_subject_2(myPaths, id)
% =========================================================================
%
% Script for basic postprocessing of the preprocessed EEG data
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================

% Load preprocessing settings
cfg = preproc_parameters;

% -------------------------------------------------------------------------
% Define paths and files
% -------------------------------------------------------------------------
subject = preproc_folders_subject(id, myPaths, 2);

% Print
fprintf('\n==================================================================\n');
fprintf('%s | %s | %s dataset | processing part 2 | pipeline v%s\n', myPaths.group, subject.id, myPaths.task, myPaths.codever);
fprintf('==================================================================\n');

% Load cleaned data
fileName = fullfile(subject.data, subject.clnfile0);
if exist(fileName, 'file') == 2
    fprintf('\n================================\n');
    fprintf('Loading data\n');
    fprintf('================================\n');

    load(fileName, 'EEG');
    fprintf('Done!\n');
else
    warning([subject.id ' is missing preprocessed ' myPaths.task ' data. Skipping...']); return;
end

% Time
t0 = datetime("now");
EEG.ALSUTRECHT.procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};

% -------------------------------------------------------------------------
% Basic post-preprocessing
% -------------------------------------------------------------------------

% Filter lowpass (should be done on continuous data)
EEG = do_filtering(EEG, 'lowpass', cfg.flt);

% % Empirical frequency boundaries (only RS)
% EEG = estimate_gedBounds(EEG);

% IAF (only RS?)
EEG = check_iaf(EEG);

% Epoch
EEGcell = epoch_data(EEG, cfg.trg);
clearvars EEG

% -------------------------------------------------------------------------
% Loop through epoched datasets
% -------------------------------------------------------------------------
for i_file = 1:length(EEGcell)
    % Extract
    EEG = EEGcell{i_file};

    % Common-average referening
    EEG = do_reref(EEG, 'aRegular');

    % Baseline correction
    EEG = do_baselinecorrection(EEG, 'traditional');

    % Detect bad epochs
    [EEG, num_trials] = detect_badepochs(EEG, cfg);

    % Final estimates
    fprintf('\n================================\n');
    fprintf('Final data estimates\n');
    fprintf('================================\n');

    % EOG leftovers
    fprintf('\nChecking EOG leftovers..\n');
    EEG = check_blink_residuals(EEG);
    num_trials(end + 1) = EEG.trials;

    % EMG leftovers
    fprintf('\nChecking EMG leftovers...\n');
    [slopes, mask_emg_matrix] = detect_emg(EEG, cfg);
    emg_leftover = mean(mask_emg_matrix, 'all');
    fprintf('EMG leftovers: %1.2f\n', emg_leftover);

    % Median voltage shift
    fprintf('\nEstimating voltage range...\n\n');
    voltage_shift = median(range(EEG.data, 2), 3);

    % Channel correlation matrix
    fprintf('Estimating channel correlation matrix..\n');
    EEG = estimate_channelcov(EEG);

    % Automagic metrics
    bad_chans = length(EEG.ALSUTRECHT.badchaninfo.badElectrodes);
    quality = report_automagic(EEG, bad_chans, 'ThresholdOHA', 30, 'ThresholdStd', 15);

    fprintf('\nAutomagic quality metrics:\n');
    fprintf('RBC (Bad Channels):   %.2f%%\n', quality.RBC * 100);
    fprintf('OHA (High Amp):       %.2f%%\n', quality.OHA * 100);
    fprintf('THV (Bad Times):      %.2f%%\n', quality.THV * 100);
    fprintf('CHV (Noisy Channels): %.2f%%\n', quality.CHV * 100);
    fprintf('Status:               %s\n',     quality.status);

    % Log
    EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal = voltage_shift;
    EEG.ALSUTRECHT.epochRejections.muscle2 = emg_leftover;
    EEG.ALSUTRECHT.leftovers.muscle2       = emg_leftover;
    EEG.ALSUTRECHT.automagicmetrics        = quality;

    % Final reports/plots
    generate_finalplots(EEG, num_trials, myPaths.task, i_file, cfg);

    % EEGLAB / BIDS metadata
    EEG = add_bidsmetadata(EEG, subject);

    % Save data
    export_data(EEG, subject, i_file);
end

% Report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf('Finished: %s\n', t1);
fprintf('Running time: %d min.\n\n', dd);

end