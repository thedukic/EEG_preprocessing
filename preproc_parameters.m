function cfg = preproc_parameters
%
% If you want to change something, discuss with the rest of the team first
%

% Preprocessing ver/run 1
cfg.rnum = '1';

%% IIR Butterworth filters
% The filters are twopass, so the filtering will be done twice
% using only the half of the defined order
%
% Input: [cutoff, final order]
% hp: highpass
% lp: lowpass

% MMN/SART (ERP) filter
cfg.flt.erp.hp = [0.25, 4];
cfg.flt.erp.lp = [80, 4];

% RS/MT filter
cfg.flt.rsmt.hp = [1, 4];
cfg.flt.rsmt.lp = [80, 4];

% EXT filter
cfg.flt.ext.hp = [0.1, 2];
cfg.flt.ext.lp = [80, 2];

% EMG filter
cfg.flt.emg.hp = [5, 4];
cfg.flt.emg.lp = [];

%% Bad channel/period detection
% Flat electrode duration
cfg.bch.flatDuration                = 4;    % [s]

% PREP: settings
cfg.bch.robustDeviationThreshold    = 5;    % default: 5
cfg.bch.highFrequencyNoiseThreshold = 5;    % default: 5
cfg.bch.correlationThreshold        = 0.4;  % default: 0.4
cfg.bch.badTimeThreshold            = 0.01; % default: 0.01

% PREP: RANSAC (computationally heavy and nondeterministic)
cfg.bch.ransacOff                   = true;
% cfg.bch.ransacCorrelationThreshold  = 0.75; % default: 0.8
% cfg.bch.ransacUnbrokenTime          = 0.5;  % default: 0.4

% % Only iff RANSAC is used
% cfg.bch.iter.num    = 10;
% cfg.bch.iter.frc    = 0.8;
% cfg.bch.iter.rejmax = 0.1;

% RELAX: Muscle activity
% Less stringent = -0.31, Middle Stringency = -0.59, More stringent = -0.72
cfg.bch.muscleSlopeThreshold        = -0.5;
cfg.bch.slopeTime                   = 0.01;

% RELAX: Slow drifts
% cfg.bch.DriftSeverityThreshold      = 12;  % default: 10 (MAD from the median of all electrodes)
% cfg.bch.driftSlopeThreshold         = -4;

% EEGLAB: ASR
cfg.bch.asr                         = 25; % recommandation: 20-30

%% ICA and ICLabels
cfg.ica.type    = 'CUDAICA';
cfg.ica.icMax   = 50; % PCA reduction prior ICA
%  {'Brain' 'Muscle' 'Eye' 'Heart' 'Line Noise' 'Channel Noise' 'Other'}
cfg.ica.iclabel = ....
    [NaN NaN; 0.8 1; 0.7 1; 0.8 1; NaN NaN; 0.7 1; NaN NaN];

%% Event triggers
cfg.trg.mmn   = {[12 17],[-0.2 0.5]};
cfg.trg.sart1 = {[3 6],[-0.2 0.9]};
cfg.trg.sart2 = {1,[-0.45 0.45]};
cfg.trg.mt    = {[21 31 51],[-5 10]};
cfg.trg.rs    = {2, 0.75}; % 2s, 0.75 overlap

end