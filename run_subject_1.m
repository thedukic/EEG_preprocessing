function run_subject_1(myPaths, id)
% RUN_SUBJECT_1 Continuous EEG preprocessing and artefact rejection (Part 1).
%
% Syntax:
%   run_subject_1(myPaths, id)
%
% Description:
%   Executes the first stage of the automated EEG preprocessing pipeline on
%   continuous 128-channel BioSemi ActiveTwo recordings.
%
%   Key processing steps:
%     1. Data import (.bdf) and event trigger alignment.
%     2. Hardware checks: CMS/DRL dropout detection and DC offset estimation.
%     3. Resampling (256 Hz), high-pass FIR filtering, and robust referencing.
%     4. Two-stage spectral line-noise suppression (50 Hz and harmonics).
%     5. External electrode derivation (bipolar VEOG, HEOG, ECG / earlobe PCA).
%     6. Spatial filtering (STAR, GEDAI) and artefact template generation.
%     7. Bad channel detection, spherical spline interpolation, and regular rereferencing.
%     8. ICA decomposition (CUDAICA / RUNICA) with subspace rank preservation.
%     9. Automated classification and rejection of ocular, cardiac, and myogenic ICs.
%    10. Export of interim preprocessed dataset and participant QA report.
%
% Inputs:
%   myPaths - Structure containing directory configurations:
%               .mycodes     : Path to pipeline source code
%               .rootrawdata : Directory containing raw BDF files
%               .rootpreproc : Output root directory for preprocessed data
%   id      - String or character array specifying participant ID (e.g., 'SUBJ001')
%
% Outputs:
%   None (processed dataset, figures, and QA logs are written directly to disk).
%
% ALS Centre, University Medical Centre Utrecht
% License: GNU General Public License v3.0

% Load preprocessing settings
cfg = preproc_parameters;

% -------------------------------------------------------------------------
% Define paths and files
% -------------------------------------------------------------------------
subject = preproc_folders_subject(id, myPaths, 1);

% Print
fprintf('==================================================================\n');
fprintf('%s | %s | %s dataset | processing part 1 | pipeline v%s\n', subject.group, subject.id, subject.task, subject.codever);
fprintf('==================================================================\n');

% -------------------------------------------------------------------------
% Basic data preparation
% -------------------------------------------------------------------------
% Find files
subject.datablocks = list_datasets(subject.rawdata, subject.task);

% Load those files
if ~isempty(subject.datablocks)
    EEG = load_biosemidata(subject, myPaths);
else
    warning([subject.id ' is missing ' subject.task ' data. Skipping...']); return;
end

% Make a folder for this participant
if exist(subject.preproc, 'dir') ~= 7
    mkdir(subject.preproc);
    mkdir(subject.data);
    mkdir(subject.figures);
    mkdir(subject.qa);
else
    warning('The folder already exists. The old and the new files might be mixed:\n%s', subject.preproc);
end

% Open a report
subject = report_proc(subject, myPaths, 'open');

% Manually fix some datasets
EEG = do_manualfix(EEG, subject, myPaths);

% Fix events
EEG = fix_events(EEG);

% Add subject/channel info
EEG = add_info(EEG, subject, cfg);

% Estimate electrode offsets
EEG = estimate_electrodeoffsets(EEG, cfg);

% Remove CMS drop-outs
EEG = detect_dropouts(EEG, cfg);

% Resample
EEG = do_resampling(EEG, 256);

% Detect flat channels
EEG = remove_flatelectrodes(EEG, cfg);

% % Remove extremely bad epochs
% [EEG, flagExclude] = remove_extremeperiods2(EEG);
% % Check if it is worth continuing
% if flagExclude, warning('%s has very noisy data. Skipping...', subject.id); return; end

% Mark where each RS block starts/ends
if ismember(upper(subject.task), {'RS', 'EO', 'EC'})
    EEG = make_blockmasks(EEG);
end

% Keep event info
EEG = extract_eventinfo(EEG, cfg.trg);

% Filter highpass only
% EEG = do_filtering_fir(EEG);
EEG = do_filtering(EEG, 'highpass', cfg.flt);

% % Electrode correlations (needs more testing)
% EEG = detect_swappedelectrodes(EEG, cfg);

% Reference
EEG = do_reref(EEG, 'aRobust');

% % Make a copy
% EEGRAW = EEG;

% -------------------------------------------------------------------------
% Reduce spectral peaks
% -------------------------------------------------------------------------
% Remove line noise
EEG = reduce_linenoise1(EEG, cfg);

% Remove line noise leftovers
EEG = reduce_linenoise2(EEG, cfg);

% % Remove other spectral peaks (needs more testing)
% EEG = reduce_spectrapeaks(EEG, cfg);

% -------------------------------------------------------------------------
% Organise
% -------------------------------------------------------------------------
% Make EXT bipolar
EEG = make_extbipolar(EEG);

% Merge blocks
EEG = merge_eeglabblocks(EEG);

% Separate EXT
[EEG, EMG, EXT] = separate_electrodetypes(EEG);

% -------------------------------------------------------------------------
% Check EC-RS blinks (needs more testing)
% -------------------------------------------------------------------------
% if ismember(upper(subject.task), {'RS', 'EC'})
%     EEG = check_eyesclosedeyeblinks(EEG, EXT, cfg);
% end

% -------------------------------------------------------------------------
% Reduce artifacts
% -------------------------------------------------------------------------
% STAR
if cfg.star.do
    EEG = do_star(EEG, 'eeg', cfg);
else
    fprintf('\nSTAR cleaning: Skipped! Turned off by the user.\n');
end

% GEDAI
if cfg.gedai.do
    EEG = do_gedai(EEG, cfg);
else
    fprintf('\nGEDAI cleaning: Skipped! Turned off by the user.\n');
end

% % CCA
% if cfg.cca.do
%     EEG = do_cca(EEG, cfg);
% else
%     fprintf('\CCA cleaning: Skipped! Turned off by the user.\n');
% end

% -------------------------------------------------------------------------
% Generate artifact templates
% -------------------------------------------------------------------------
generate_ictemplateweights(EEG, EMG, EXT, cfg);

% -------------------------------------------------------------------------
% Remove bad EEG channels
% -------------------------------------------------------------------------
% Remove noisy electrodes
EEG = remove_noisyelectrodes(EEG, cfg);

% Report bad/removed electrodes
EEG = report_badelectrodes(EEG, cfg);

% -------------------------------------------------------------------------
% Deeper cleaning
% -------------------------------------------------------------------------
% % MWF
% if cfg.mwf.do
%     EEG = do_mwf(EEG,EXT,cfg);
% else
%     fprintf('\nMWF cleaning: Skipped! Turned off by the user.\n');
% end
% % ASR
% if cfg.asr.do
%     EEG = do_asr(EEG,cfg);
% else
%     fprintf('\nASR cleaning: Skipped! Turned off by the user.\n');
% end

% Interpolate bad electrodes
EEG = do_channelinterp(EEG, 'spherical');

% Reference
EEG = do_reref(EEG, 'aRegular');

% ICA
% EEG = do_relica(EEG, subject);
EEG = do_ica(EEG, cfg);

% Detect artifact ICs
% EEG = detect_ic_bad(EEG, EXT, EMG, cfg); % not working
EEG = detect_badcomponents(EEG, EXT, EMG, cfg);

% Remove artifact ICs
EEG = remove_badcomponents(EEG, cfg);

% Report ICA
EEG = report_ica(EEG, cfg);

% -------------------------------------------------------------------------
% Clean EMG (not done here?)
% -------------------------------------------------------------------------
% if ismember(upper(subject.task), 'MT')
%     EMG = clean_emg(EMG);
%     EMG2 = do_cca(EMG, cfg);
%
%     % EMG.etc = rmfield(EMG.etc, 'clean_sample_mask');
%     % EMG2.etc = rmfield(EMG2.etc, 'clean_sample_mask');
%     % vis_artifacts(EMG2, EMG);
% end

% -------------------------------------------------------------------------
% Finalise
% -------------------------------------------------------------------------
% Merge channels/sets
EEG = merge_electrodetypes(EEG, EMG, EXT);
clearvars EXT EMG

% Report artifact leftovers
EEG = report_leftovers(EEG, 1, cfg);

% Interim data saving
% EEG = pop_saveset(EEG, 'filename', [subject.id '_' subject.visit '_' subject.task '_tmp.set'], 'filepath', subject.preproc);

% % Visually check the cleaning
% compare_visually(EEG, EEGRAW, cfg.trg);

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
% Close the report
subject = report_proc(subject, myPaths, 'close');

% Save data
export_data(EEG, subject, 0);

end