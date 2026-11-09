function run_subject_1(myPaths, id)
% RUN_SUBJECT_1 Continuous EEG preprocessing and artefact rejection (Part 1).
%
% Syntax:
%   run_subject_1(myPaths, id)
%
% Description:
%   Executes the first stage of the automated EEG preprocessing pipeline on
%   continuous EEG BioSemi ActiveTwo recordings.
%
% Inputs:
%   myPaths - Structure containing directory configurations
%   id      - String or character array specifying participant ID
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
subject = preproc_folders_subject(id, myPaths);

% Print
fprintf('==================================================================\n');
fprintf('%s | %s | %s dataset | processing part 1 | pipeline v%s\n', subject.group, subject.id, subject.task, subject.codever);
fprintf('==================================================================\n');

% -------------------------------------------------------------------------
% Basic data preparation
% -------------------------------------------------------------------------
% Find and load files
EEG = load_biosemidata(subject, cfg);

if isempty(EEG)
    return;
else
    clearvars subject
end

% Open a report
EEG = report_proc(EEG, myPaths, 'open');

% Manually fix some datasets
EEG = do_manualfix(EEG);

% Fix events
EEG = fix_events(EEG);

% Remove CMS drop-outs
EEG = detect_dropouts(EEG, cfg);

% Remove extremely bad segments
% [EEG, flag_exclude] = remove_extremeperiods2(EEG);
[EEG, flag_exclude] = detect_extremes(EEG);

% Check if it is worth continuing
if flag_exclude, warning('%s has very noisy data. Skipping...', subject.id); return; end

% Resample
EEG = do_resampling(EEG, 256);

% Detect flat channels
EEG = remove_flatelectrodes(EEG, cfg);

% Mark where each RS block starts/ends
if ismember(upper(EEG(1).ALSUTRECHT.subject.task), {'RS', 'EO', 'EC'})
    EEG = make_blockmasks(EEG);
end

% Keep event info
EEG = extract_eventinfo(EEG, cfg);

% Filter highpass only
EEG = do_filtering(EEG, 'highpass', 'eeg', cfg);

% % Electrode correlations (needs more testing)
% EEG = detect_swappedelectrodes(EEG, cfg);

% Reference
EEG = do_reref(EEG, 'aRobust');

% % Make a copy
% EEGRAW = EEG;

% Save raw power for later
save_raw_power(EEG, cfg);

% -------------------------------------------------------------------------
% Reduce line noise
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
% Sparse artefact reduction
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

% % ASR
% if cfg.asr.do
%     EEG = do_asr(EEG,cfg);
% else
%     fprintf('\nASR cleaning: Skipped! Turned off by the user.\n');
% end

% % MWF
% if cfg.mwf.do
%     EEG = do_mwf(EEG,EXT,cfg);
% else
%     fprintf('\nMWF cleaning: Skipped! Turned off by the user.\n');
% end

% -------------------------------------------------------------------------
% Bad Channel Rejection
% -------------------------------------------------------------------------
% Detect and remove
EEG = remove_noisyelectrodes(EEG, cfg);

% Report
report_badelectrodes(EEG, cfg);

% -------------------------------------------------------------------------
% Generate artifact templates
% -------------------------------------------------------------------------
% Reference
EEG = do_reref(EEG, 'aRegular');

% Make templates
generate_ictemplateweights(EEG, EXT, cfg);

% -------------------------------------------------------------------------
% ICA
% -------------------------------------------------------------------------
% % Interpolate bad electrodes
% EEG = do_channelinterp(EEG, 'spherical');
%
% % Reference
% EEG = do_reref(EEG, 'aRegular');

% ICA
% EEG = do_relica(EEG, subject);
EEG = do_ica(EEG, cfg);

% Interpolate bad electrodes
EEG = do_channelinterp(EEG, 'spherical');

% ICA: Detect artifacts
EEG = detect_badcomponents_new(EEG, EXT, EMG, cfg);

% ICA: Remove artifacts
EEG = remove_badcomponents(EEG, cfg);

% ICA: Report
report_ica(EEG, cfg);

% -------------------------------------------------------------------------
% Deeper cleaning
% -------------------------------------------------------------------------
% GBF
if cfg.gbf.do
    EEG = do_gbf_time(EEG, cfg);
else
    fprintf('\GBF cleaning: Skipped! Turned off by the user.\n');
end

% % CCA
% if cfg.cca.do
%     EEG = do_cca(EEG, cfg);
% else
%     fprintf('\CCA cleaning: Skipped! Turned off by the user.\n');
% end

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

% % MWF / DSS
% if DATA.ALSUTRECHT.leftovers.blink1.flag_redo
%     EEG = do_mwf_blink(EEG);
% end

% Interim data saving
% EEG = pop_saveset(EEG, 'filename', [subject.id '_' subject.visit '_' subject.task '_tmp.set'], 'filepath', subject.preproc);

% % Visually check the cleaning
% compare_visually(EEG, EEGRAW, cfg.trg);

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
% Close the report
EEG = report_proc(EEG, myPaths, 'close');

% Save data
save_data(EEG, 0);

end