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
cfg.flt.erp.hp = [0.5, 4];
cfg.flt.erp.lp = [80, 4];

% MT filter
cfg.flt.mt.hp = [1, 4];
cfg.flt.mt.lp = [80, 4];

% RS filter
cfg.flt.rs.hp = [0.5, 4]; % 1?
cfg.flt.rs.lp = [80, 4];

% EXT filter
cfg.flt.ext.hp = [0.15, 4];
cfg.flt.ext.lp = [80, 4];

% EMG filter
cfg.flt.emg.hp = [5, 4];
cfg.flt.emg.lp = [];

%% Bad channel/period detection
% Flat electrode duration
cfg.bch.flatDuration                = 4;    % default: 4 [s]

% PREP: settings
cfg.bch.robustDeviationThreshold    = 5;    % default: 5
cfg.bch.highFrequencyNoiseThreshold = 5;    % default: 5
cfg.bch.correlationThreshold        = 0.40; % default: 0.4
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
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7590828
% If 7-45 then use muscleSlopeThreshold = 0 ?
% If 7-70 the use the above muscleSlopeThreshold, like -0.31 ?
cfg.bch.muscleSlopeThreshold        = 0;
cfg.bch.muscleSlopeFreq             = [7 70]; % 7-45 or 7-70/75
cfg.bch.muscleSlopeTime             = 0.50;
cfg.bch.maxProportionOfBadElec      = 0.15;  % that can be deleted using this EMG detection

% RELAX: Slow drifts
% cfg.bch.DriftSeverityThreshold      = 12;  % default: 10 (MAD from the median of all electrodes)
% cfg.bch.driftSlopeThreshold         = -4;

% EEGLAB: ASR
% cfg.bch.asr                         = 20;    % recommandation: 20-30

%% ICA and ICLabels

% ICA algorithm
gpuCount = gpuDeviceCount;
if gpuCount > 0
    % disp('GPU is available. The preprocessing will be quicker as CUDAICA will be used.');
    cfg.ica.type1    = 'CUDAICA';
    cfg.ica.type2    = 'CUDAICA'; % AMICA
else
    % disp('No GPU available. The preprocessing will be slower as RUNICA will be used instead of CUDAICA.');
    cfg.ica.type1    = 'RUNICA';
    cfg.ica.type2    = 'RUNICA'; % AMICA
end

% PCA reduction prior ICA
cfg.ica.icMax   = 70;

% ICLabel likelihoods
cfg.ica.iclabel = ...
    [NaN NaN;  % Brain
    0.4 1;     % Muscle - used
    0.4 1;     % Eye (VEOG+HEOG) - used
    NaN NaN;   % Heart - not used
    NaN NaN;   % Line noise - not needed
    0.4 1;     % Channel noise - used
    NaN NaN];  % Other

% cfg.ica.blinkchans = {'C8','C9','C10','C14','C15','C16','C17','C18','C19','C27','C28','C29','C30','C31','C32','C26','C20','C13','C21'};
cfg.ica.blinkchans = {'C8','C9','C10','C14','C15','C16','C17','C18','C19','C27','C28','C29','C30','C31','C32'};
% BlinkElectrodes = {'C14','C15','C16','C17','C18','C19','C27','C28','C29'};
% BlinkElectrodes = {'C8','C17','C29','C30'};

%% Event triggers
cfg.trg.mmn   = {[12 17],[-0.2 0.5]};
cfg.trg.sart1 = {[3 6],[-0.2 0.9]};
cfg.trg.sart2 = {[1 11],[-0.4 0.4]};
cfg.trg.mt    = {[21 31 51],[-5 10]};
cfg.trg.rs    = {2,0.75};             % 2s, 0.75 overlap

%% Bad epoch rejection
cfg.epoch.rejectAmp                            = 75;  % even >75 is too low for ALS30344
% cfg.epoch.singleChannelImprobableDataThreshold = 5; % MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
% cfg.epoch.allChannelImprobableDataThreshold    = 3; % SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
% cfg.epoch.singleChannelKurtosisThreshold       = 5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis
% cfg.epoch.allChannelKurtosisThreshold          = 3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis
cfg.epoch.singleChannelImprobableDataThreshold = 8; % MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
cfg.epoch.allChannelImprobableDataThreshold    = 6; % SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
cfg.epoch.singleChannelKurtosisThreshold       = 8; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis
cfg.epoch.allChannelKurtosisThreshold          = 6; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis

end