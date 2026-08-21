function cfg = preproc_parameters
% =========================================================================
% PREPROCESSING PIPELINE PARAMETER CONFIGURATION
% ALS Centre, University Medical Centre Utrecht
% =========================================================================
%
% Script for defining global preprocessing parameters across experimental
% paradigms (Resting-State, MMN, SART, and Motor Task).
%
% -------------------------------------------------------------------------
% IMPORTANT NOTICE
% -------------------------------------------------------------------------
% Parameter modifications alter cohort-level reproducibility. Any changes
% must be discussed and agreed upon with the rest of the team first.
% =========================================================================


% =========================================================================
% BIOSEMI 128 (ABCD) REGION OF INTEREST (ROI) CONFIGURATIONS
% =========================================================================

% 1. Cognitive ERP Hubs
% -------------------------------------------------------------------------
% MMN: Frontocentral negativity (Centred around Fz / FCz)
% Midline Fz (C21) + surrounding anterior/posterior and bilateral C/D leads
cfg.roi.mmn = { ...
    'C21', ...              % Fz (Midline core)
    'C20', 'C19', ...       % Anterior to Fz
    'C22', 'C23', ...       % Posterior to Fz (towards FCz / Cz)
    'C25', 'C26', ...       % Left frontocentral
    'C12', 'C13', ...       % Right frontocentral
    'D1', 'D2', ...         % Left vertex-adjacent ring
    'C1', 'C2'              % Right vertex-adjacent ring
    };

% SART: Centroparietal P300 / P3b (Centred around Pz)
% Midline Pz (A19) + surrounding centroparietal A-leads
cfg.roi.sart_p300 = { ...
    'A19', ...              % Pz (Midline core)
    'A4', 'A3', 'A2', ...   % Anterior to Pz (towards CPz / Cz)
    'A20', 'A21', ...       % Posterior to Pz (towards POz)
    'A5', 'A18', ...        % Left centroparietal
    'A32', 'A31'            % Right centroparietal
    };


% 2. Spectral & Sensorimotor Hubs
% -------------------------------------------------------------------------
% Resting-State: Occipital Alpha Hub (Centred around Oz)
% Midline Oz (A23) + surrounding occipital / parieto-occipital ring
cfg.roi.rs_alpha = { ...
    'A23', ...              % Oz (Occipital core)
    'A22', 'A21', ...       % Anterior to Oz (POz)
    'A24', 'A25', ...       % Inferior / Inion
    'A15', 'A16', 'A17', ...% Left occipital / O1
    'A28', 'A29', 'A30', ...% Right occipital / O2
    'A14', 'A27'            % Outer inferolateral occipital
    };

% Motor Task: Left Sensorimotor Core (Centred around C3 for Right-Hand MT)
% C3 (D19) + immediate surrounding ring in D-Bank
cfg.roi.motor_left = { ...
    'D19', ...              % C3 (Core)
    'D18', 'D14', ...       % Medial to C3
    'D20', 'D21', ...       % Lateral to C3
    'D12', 'D11', 'D13', ...% Anterior to C3 (Premotor / FC3)
    'D17', 'D28', 'D27'     % Posterior to C3 (Sensory / CP3)
    };

% Motor Task: Right Sensorimotor Core (Centred around C4 for Left-Hand MT)
% C4 (B22) + immediate surrounding ring in B-Bank
cfg.roi.motor_right = { ...
    'B22', ...              % C4 (Core)
    'B21', 'B20', ...       % Medial to C4
    'B23', 'B24', ...       % Lateral to C4
    'B31', 'B30', 'B32', ...% Anterior to C4 (Premotor / FC4)
    'B19', 'B18', 'B17'     % Posterior to C4 (Sensory / CP4)
    };

% Combined bilateral mask for Motor task validation
cfg.roi.motor_all = [cfg.roi.motor_left, cfg.roi.motor_right];

% 3. Artifact Validation & Template Masks
% -------------------------------------------------------------------------
% Vertical EOG / Blink ROI (Anterior Pole surrounding Fpz / C17)
cfg.roi.blink = { ...
    'C17', ...              % Fpz (Anterior pole apex)
    'C18', 'C19', ...       % Midline supra-orbital
    'C29', 'C30', ...       % Left prefrontal outermost (Circle 8)
    'C16', 'C8', ...        % Right prefrontal outermost (Circle 8)
    'C28', 'C27', ...       % Left prefrontal inner (Circles 7, 6)
    'C15', 'C14'            % Right prefrontal inner (Circles 7, 6)
    };

% Horizontal EOG / Saccade: Left Outer Fronto-Temporal Rim (F7 equivalent)
cfg.roi.saccade_left = { ...
    'D7', ...               % Core Left Fronto-Lateral Peak (F7 / LO1)
    'D6', 'D8', ...         % Flanking leads
    'C30', 'C31', ...       % Anterior orbital border
    'D5', 'D9'              % Inner fronto-temporal perimeter
    };

% Horizontal EOG / Saccade: Right Outer Fronto-Temporal Rim (F8 equivalent)
cfg.roi.saccade_right = { ...
    'C7', ...               % Core Right Fronto-Lateral Peak (F8 / LO2)
    'C6', 'B27', ...        % Flanking leads
    'C8', 'C9', ...         % Anterior orbital border
    'C5', 'B28'             % Inner fronto-temporal perimeter
    };

% Combined bilateral mask for HEOG dipole validation
cfg.roi.saccade_all = [cfg.roi.saccade_left, cfg.roi.saccade_right];

% =========================================================================
% IIR Butterworth Filter Configuration
% =========================================================================
% Zero-phase application via filtfilt (effective order is doubled).
%
% Format : [cutoff_hz, effective_order]
% hp     : High-pass
% lp     : Low-pass ([] denotes open / unconstrained)

% -------------------------------------------------------------------------
% 1. Cognitive ERP / ERSP (MMN / SART)
% -------------------------------------------------------------------------
% 0.3 Hz preserves P300/MMN and theta without baseline drift.
% 60 Hz (effective order 8) preserves flat response to 50 Hz while cutting >75 Hz EMG.
cfg.flt.erp.hp = [0.3, 4];
cfg.flt.erp.lp = [60, 8];

% -------------------------------------------------------------------------
% 2. Resting-State (RS)
% -------------------------------------------------------------------------
% 0.5 Hz removes galvanic sweat drift.
% Open LP ([]) prevents filter roll-off curvature from biasing 1/f aperiodic fitting.
cfg.flt.rs.hp = [0.5, 4];
cfg.flt.rs.lp = [];

% -------------------------------------------------------------------------
% 3. Motor Task EEG (ERSP / CMC)
% -------------------------------------------------------------------------
% 1.0 Hz removes movement and mechanical cable sway.
% 60 Hz preserves the 15-45 Hz CMC band with <0.4 dB attenuation at 45 Hz.
cfg.flt.mt.hp = [1.0, 4];
cfg.flt.mt.lp = [60, 8];

% -------------------------------------------------------------------------
% 4. Motor Task EMG (CMC)
% -------------------------------------------------------------------------
% 10 Hz preserves descending 15-30 Hz raw beta drive without motion artifacts.
% Open LP ([]) retains broad motor unit action potential profiles.
cfg.flt.emg.hp = [10, 4];
cfg.flt.emg.lp = [];

% -------------------------------------------------------------------------
% 5. External Channels
% -------------------------------------------------------------------------
% % Fitler settings are matched inside the pipeline to EEG to maintain phase alignment (ICA).
% cfg.flt.ext.hp = [0.3, 4];
% cfg.flt.ext.lp = [60, 8];

% =========================================================================
% Bad channel/data detection
% =========================================================================
% cfg.emg.asr_z_low = -3.5;
% cfg.emg.asr_z_high = 7.0;

cfg.emg.slope_time        = 0.50;
cfg.emg.slope_freq_1      = [7 75]; % [7 45] / [7 75]
cfg.emg.slope_threshold_1 = 0;
cfg.emg.slope_freq_2      = [40 70];
cfg.emg.slope_threshold_2 = -0.2;
cfg.emg.slope_freq_excl   = 50;

% =========================================================================
% Bad channel/data detection
% =========================================================================
% Flat electrode duration
cfg.channel.flatDuration = 4;    % default: 4 [s]

% PREP: Settings
% badTimeThreshold = 0.01 (PREP default) is overly aggressive: it declares a channel
% bad if it drops below threshold for just 1% of time (7 s in a 12 min run), causing
% healthy channels with brief swallows/muscle bursts to be over-interpolated.
%
% R = 0.40 is appropriate for high-density caps due to dense electrode spacing (~2 cm).
% If peripheral/temporal rim channels are falsely flagged due to one-sided neighbour interpolation, drop to 0.35.
cfg.channel.robustDeviationThreshold    = 5;    % default: 5
cfg.channel.highFrequencyNoiseThreshold = 5;    % default: 5
cfg.channel.correlationThreshold        = 0.40; % default: 0.4; higher values -> more stringent
cfg.channel.badTimeThreshold            = 0.05; % default: 0.01

% PREP: RANSAC (computationally heavy and nondeterministic)
cfg.channel.ransacOff                   = true;
% cfg.channel.ransacCorrelationThreshold  = 0.75; % default: 0.8
% cfg.channel.ransacUnbrokenTime          = 0.5;  % default: 0.4
% cfg.channel.iter.num    = 10;
% cfg.channel.iter.frc    = 0.8;
% cfg.channel.iter.rejmax = 0.1;

% RELAX: Muscle activity
% If 7-45 Hz, -> use emgSlopeThreshold = 0
% If 7-70 Hz, -> Less stringent = -0.31, Middle Stringency = -0.59, More stringent = -0.72
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7590828
% https://mne.tools/stable/generated/mne.preprocessing.ICA.html
% cfg.channel.emgSlopeFreq           = [7 45];
% cfg.channel.emgSlopeThreshold      = 0;
% cfg.channel.emgSlopeTime           = 0.50;
cfg.channel.prop_badchan_max = 0.10;   % 10% (~128/10)of total can be deleted using this EMG detection only

% =========================================================================
% Cleaning methods
% =========================================================================
% STAR
cfg.star.do = true;

% GEDAI
cfg.gedai.do = true;

% CCA
cfg.cca.do = false;

% MWF
cfg.mwf.do = false;

% ASR
cfg.asr.do = false;
% cfg.asr.std = 20; % recommandation: 20-30

% =========================================================================
% ICA
% =========================================================================
% ICA algorithm
gpuCount = gpuDeviceCount;
if gpuCount > 0
    % disp('GPU is available. The preprocessing will be quicker as CUDAICA will be used.');
    cfg.ica.type = 'CUDAICA'; % AMICA
else
    % disp('No GPU available. The preprocessing will be slower as RUNICA will be used instead of CUDAICA.');
    cfg.ica.type = 'RUNICA'; % AMICA
end

% PCA reduction prior to ICA
% cfg.ica.num_ica = [70 50]; %  (old setting)
cfg.ica.num_ica = [100 80 50];

% % ICLabel (old)
% cfg.ica.iclabel = [ ...
%     NaN NaN;   % Brain
%     0.6 1;     % Muscle
%     0.6 1;     % Eye (VEOG + HEOG)
%     0.6 1;     % Heart - not very good
%     NaN NaN;   % Line noise - not needed
%     0.6 1;     % Channel noise
%     NaN NaN];  % Other

% Base threshold template for all paradigms
cfg.ica.iclabel.base = [ ...
    NaN  NaN;    % 1. Brain
    0.80 1.0;    % 2. Muscle (placeholder)
    0.70 1.0;    % 3. Eye
    0.50 1.0;    % 4. Heart
    NaN  NaN;    % 5. Line noise
    0.80 1.0;    % 6. Channel noise
    NaN  NaN];   % 7. Other

% Muscle thresholds per paradigm
cfg.ica.iclabel.muscle_thresh.erp = 0.85; % MMN, SART
cfg.ica.iclabel.muscle_thresh.mt  = 0.85; % Motor task
cfg.ica.iclabel.muscle_thresh.rs  = 0.65; % Resting-State

% icablinkmetrics
% cfg.ica.blinkchans = {'C8','C9','C10','C14','C15','C16','C17','C18','C19','C27','C28','C29','C30','C31','C32','C26','C20','C13','C21'};
% cfg.ica.blinkchans = {'C8','C9','C10','C14','C15','C16','C17','C18','C19','C27','C28','C29','C30','C31','C32'};
% cfg.ica.blinkchans = {'C14','C15','C16','C17','C18','C19','C27','C28','C29'};
% cfg.ica.blinkchans = {'C8','C17','C29','C30'};
% cfg.ica.blinkchans = {'C29','C17','C16'};
cfg.ica.blinkchans = cfg.roi.blink;

% Muscle ICs
cfg.ica.emg = true;
% cfg.ica.emgSlopeFreq1      = [7 45];
% cfg.ica.emgSlopeThreshold1 = 0;
% cfg.ica.emgSlopeFreq2      = [40 70];
% cfg.ica.emgSlopeThreshold2 = -0.2;

% % Cutoff: ICs will be filtered out above this freq (i.e., keeps low freqs)
% % https://www.biorxiv.org/content/10.1101/2024.06.06.597688v1.full.pdf
% cfg.ica.emgfilter = 15; % Suggested: 15 Hz

% Channel ICs
cfg.ica.channel = true;

% Generally bad components [without a precise label]
cfg.ica.bad = true;

% =========================================================================
% Event Triggers and Epoch Settings
% Format for Task:           { [trigger_ids], [t_start t_stop] (s) }
% Format for Resting-state:  { win_length (s), overlap_ratio }
% =========================================================================

% --- Event-related potentials & task paradigms ---
cfg.trg.mmn   = {[12 17],    [-0.2  0.5]};   % Mismatch Negativity
cfg.trg.sart1 = {[3 6],      [-1.1  1.1]};   % SART 1 (old: [-0.2 0.9])
cfg.trg.sart2 = {1,          [-1.2  1.0]};   % SART 2
cfg.trg.mt    = {[21 31 51], [-5.0 10.0]};   % Motor task

% --- Resting-state continuous segmentation ---
cfg.trg.rs1   = {2, 0.5};                    % 2 s windows, 50% overlap
% cfg.trg.rs2 = {5, 0.8};                    % 5 s windows, 80% overlap

% =========================================================================
% Bad epoch rejection
% =========================================================================
% Maximum amplitude
cfg.epoch.amplitude_max = 75;

% Interpolate EMG-contaminated trials
cfg.epoch.interpolation = 'skip'; % yes / no / skip
% If YES, then define max. number of contaminated channels:
% -> if more than this number, then the trial is removed
% -> otherwise, it is interpolated
cfg.epoch.interpolation_max = 3;

% EEGLAB rejections
% 1. Relax single-channel bounds so 1 or 2 peripheral electrodes cannot kill the epoch
cfg.epoch.singleChannelImprobableDataThreshold = 6.0; % MAD
cfg.epoch.singleChannelKurtosisThreshold       = 6.0; % SD

% 2. Keep all-channel bounds generous to only catch massive whole-head artifacts
cfg.epoch.allChannelImprobableDataThreshold    = 4; % SD across all channels
cfg.epoch.allChannelKurtosisThreshold          = 4; % SD across all channels

% % EEGLAB rejections (RELAX toolbox settings)
% cfg.epoch.singleChannelImprobableDataThreshold = 5; % MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
% cfg.epoch.allChannelImprobableDataThreshold    = 3; % SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
% cfg.epoch.singleChannelKurtosisThreshold       = 5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis
% cfg.epoch.allChannelKurtosisThreshold          = 3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis

% =========================================================================
% Figure settings
% =========================================================================
cfg.figure.visible = 'on';

end