function simulate_eeg(path_leadfield, path_save, id)
% path_leadfield = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\individual\leadfields_new\ALS02347_1_leadfield_all.mat';
% path_leadfield = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\individual\leadfields_new\ALS08518_1_leadfield_all.mat';
% path_leadfield = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\individual\leadfields_new\ALS08549_1_leadfield_all.mat';

path_save = fullfile(path_save, id);
if ~exist(path_save, 'dir'); mkdir(path_save); end

% Simulation timing
fs = 512;
duration = 120; % 2 minutes
N_samples = fs * duration;
t = (0:N_samples-1) / fs;
file_conditions = {'EO', 'EO', 'EO', 'EC', 'EC', 'EC'};

% Load FieldTrip leadfield
load(path_leadfield, 'leadfield');
valid_dipoles = find(leadfield.inside);
N_sources = length(valid_dipoles);
N_eeg = 128;

% Extract leadfield matrix L (128 x N_sources)
L = zeros(N_eeg, N_sources);
for i = 1:N_sources
    L(:, i) = leadfield.leadfield{valid_dipoles(i)}(:, 1);
end

% Fast low-rank projection matrix for cortical 1/f background
C_scalp = L * L';
[V, D] = eig(C_scalp);
D = max(real(diag(D)), 0);
A_proj = V * diag(sqrt(D)); % 128 x 128

% Extract leadfield vectors for focal oscillations
occipital_idx = round(N_sources * 0.8);
motor_idx     = round(N_sources * 0.5);
L_alpha = L(:, occipital_idx);
L_beta  = L(:, motor_idx);

% Channel labels (128 scalp + 8 EXG + 1 Status = 137 total)
labels_scalp = cell(N_eeg, 1);
prefixes = {'A', 'B', 'C', 'D'};
idx_lbl = 1;
for p = 1:length(prefixes)
    for n = 1:32
        labels_scalp{idx_lbl} = sprintf('%s%d', prefixes{p}, n);
        idx_lbl = idx_lbl + 1;
    end
end
% labels_exg = {'ECG', 'ECG R', 'VEOG S', 'VEOG I', 'HEOG L', 'HEOG R', 'Mastoid L', 'Mastoid R'}';
labels_exg = {'EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8'}';
all_labels = [labels_scalp; labels_exg; {'Status'}];
N_total    = length(all_labels);

% 1. Load and Organise ICA Topography Templates (128 x 1 each)
load('wBlink.mat',   'Blinkweights');
load('wSaccade.mat', 'Saccadeweights');
load('wHeart.mat',   'Heartweights');

topo_blink   = Blinkweights(1:128, 1);
topo_saccade = Saccadeweights(1:128, 1);
topo_heart   = Heartweights(1:128, 1);

topo_blink   = topo_blink   / max(abs(topo_blink));
topo_saccade = topo_saccade / max(abs(topo_saccade));
topo_heart   = topo_heart   / max(abs(topo_heart));

% 2. Generate 1D Artefact Activation Time Courses
heart_rate_bpm = 68;
ecg_activation = generate_synthetic_ecg(N_samples, fs, heart_rate_bpm);

blink_activation = zeros(1, N_samples);
blink_times = round(fs * 2.5 : fs * 3.8 : N_samples - fs * 2);
blink_len   = round(0.35 * fs);
blink_shape = hanning(blink_len)' * 120;
for b = 1:length(blink_times)
    idx = blink_times(b) : blink_times(b) + blink_len - 1;
    blink_activation(idx) = blink_activation(idx) + blink_shape;
end

saccade_activation = zeros(1, N_samples);
saccade_times = round(fs * 4 : fs * 6 : N_samples - fs * 4);
for s = 1:length(saccade_times)
    step_amp = (rand() - 0.5) * 60;
    idx = saccade_times(s) : min(saccade_times(s) + round(fs * 1.8), N_samples);
    saccade_activation(idx) = saccade_activation(idx) + step_amp;
end
saccade_activation = smoothdata(saccade_activation, 'gaussian', round(0.04 * fs));

[b_emg, a_emg] = butter(2, [30 150] / (fs / 2), 'bandpass');
emg_noise      = filtfilt(b_emg, a_emg, randn(1, N_samples)) * 5;
emg_envelope   = smooth(randn(1, N_samples), fs)' > 0.45;
emg_active     = emg_noise .* emg_envelope;
temporal_chans = [100, 101, 110, 111];

% 3. Main Simulation and Biosig Export Loop
for file_idx = 1:length(file_conditions)
    condition = file_conditions{file_idx};
    run_idx = sum(strcmp(file_conditions(1:file_idx), condition));

    % Fast 128-channel sensor-space 1/f background
    pink_128 = generate_pink_noise(N_samples, N_eeg)' * 2;
    data_scalp = A_proj * pink_128;

    % Occipital Alpha (10 Hz)
    if strcmp(condition, 'EC')
        alpha_amp = 15;
    else
        alpha_amp = 3;
    end
    alpha_wave = alpha_amp * sin(2 * pi * 10 * t) .* (1 + 0.2 * randn(1, N_samples));
    data_scalp = data_scalp + (L_alpha * alpha_wave);

    % Motor Beta (20 Hz)
    beta_amp  = 5;
    beta_wave = beta_amp * sin(2 * pi * 20 * t) .* (1 + 0.1 * randn(1, N_samples));
    data_scalp = data_scalp + (L_beta * beta_wave);

    % Project ICA artefacts onto scalp channels
    data_scalp = data_scalp + (topo_blink   * blink_activation);
    data_scalp = data_scalp + (topo_saccade * saccade_activation);
    data_scalp = data_scalp + (topo_heart   * (ecg_activation * 0.08));

    % Add temporal EMG bursts
    for ch = temporal_chans
        data_scalp(ch, :) = data_scalp(ch, :) + emg_active;
    end

    % Synthesise 8 external channels
    mastoid_l = data_scalp(100, :) * 0.4 + generate_pink_noise(N_samples, 1)' * 2 + ecg_activation * 0.02;
    mastoid_r = data_scalp(110, :) * 0.4 + generate_pink_noise(N_samples, 1)' * 2 + ecg_activation * 0.02;

    data_exg = [
        ecg_activation                          + generate_pink_noise(N_samples, 1)' * 8;  % ECG
        -0.35 * ecg_activation                   + generate_pink_noise(N_samples, 1)' * 4;  % ECG R
        blink_activation                        + generate_pink_noise(N_samples, 1)' * 3;  % VEOG S
        -0.3 * blink_activation                  + generate_pink_noise(N_samples, 1)' * 3;  % VEOG I
        saccade_activation                      + generate_pink_noise(N_samples, 1)' * 2;  % HEOG L
        -saccade_activation                      + generate_pink_noise(N_samples, 1)' * 2;  % HEOG R
        mastoid_l;                                                                         % Mastoid L
        mastoid_r                                                                          % Mastoid R
        ];

    % Channel 137: BioSemi Status Channel with discrete triggers
    status_channel = zeros(1, N_samples);
    if strcmp(condition, 'EO')
        trig_val = 1;
    else
        trig_val = 2;
    end
    trig_pos  = round(1.0 * fs);
    pulse_len = round(0.05 * fs);
    status_channel(trig_pos : trig_pos + pulse_len - 1) = trig_val;

    % Complete 137-channel matrix
    data_full = [data_scalp; data_exg; status_channel];

    % Export to 24-bit BDF via Biosig
    file_save = sprintf('%s_T1_%s%d.bdf', id, condition, run_idx);
    filepath = fullfile(path_save, file_save);

    record_dur = 1;
    HDR = [];
    HDR.TYPE       = 'BDF';
    HDR.FileName   = filepath;
    HDR.SampleRate = fs;
    HDR.NS         = N_total;
    HDR.Label      = all_labels;
    HDR.Dur        = record_dur;
    HDR.SPR        = fs * record_dur;
    HDR.AS.SPR     = repmat(HDR.SPR, 1, N_total);
    HDR.NRec       = floor(N_samples / HDR.SPR);
    HDR.T0         = clock;

    % Microvolt range for analogue channels (1 to 136)
    HDR.PhysMin     = repmat(-262144, 1, N_total);
    HDR.PhysMax     = repmat( 262143, 1, N_total);
    HDR.DigMin      = repmat(-8388608, 1, N_total);
    HDR.DigMax      = repmat( 8388607, 1, N_total);
    HDR.PhysDim     = repmat({'uV'}, 1, N_total);
    HDR.PhysDimCode = repmat(4275, 1, N_total);

    % Configure Status channel (Channel 137): exact 16-bit integer range
    idx_status = N_total;
    HDR.PhysMin(idx_status)     = 0;
    HDR.PhysMax(idx_status)     = 65535;
    HDR.DigMin(idx_status)      = 0;
    HDR.DigMax(idx_status)      = 65535;
    HDR.PhysDim{idx_status}     = 'count';
    HDR.PhysDimCode(idx_status) = 512; % IEEE 11073 code for dimensionless/count

    HDR = sopen(HDR, 'w');
    total_samples = HDR.NRec * HDR.SPR;
    HDR = swrite(HDR, double(data_full(:, 1:total_samples)'));
    HDR = sclose(HDR);

    fprintf('Exported %s with 137 channels (128 scalp + 8 EXG + Status)\n', filepath);
end
end

% Helper: Synthetic P-Q-R-S-T ECG Generator
function ecg = generate_synthetic_ecg(N_samples, fs, bpm)
ecg = zeros(1, N_samples);
rr_interval = round(fs * (60 / bpm));

params.P = [  45, -0.160, 0.035];
params.Q = [ -60, -0.035, 0.012];
params.R = [ 850,  0.000, 0.016];
params.S = [-180,  0.040, 0.016];
params.T = [ 190,  0.240, 0.065];

beat_centers = round(fs * 0.6) : rr_interval : (N_samples - fs * 0.6);
win_samples  = round(0.4 * fs);

for b = 1:length(beat_centers)
    c = beat_centers(b);
    idx = max(1, c - win_samples) : min(N_samples, c + win_samples);
    t_rel = (idx - c) / fs;

    beat_wave = params.P(1) * exp(-((t_rel - params.P(2)).^2) / (2 * params.P(3)^2)) + ...
        params.Q(1) * exp(-((t_rel - params.Q(2)).^2) / (2 * params.Q(3)^2)) + ...
        params.R(1) * exp(-((t_rel - params.R(2)).^2) / (2 * params.R(3)^2)) + ...
        params.S(1) * exp(-((t_rel - params.S(2)).^2) / (2 * params.S(3)^2)) + ...
        params.T(1) * exp(-((t_rel - params.T(2)).^2) / (2 * params.T(3)^2));

    ecg(idx) = ecg(idx) + beat_wave;
end
end

% Helper: Pink Noise Generator
function y = generate_pink_noise(N_samples, N_channels)
X = fft(randn(N_samples, N_channels));
num_unique_pts = floor(N_samples / 2) + 1;
f = (0:num_unique_pts - 1)';
amp_scale = [0; 1 ./ sqrt(f(2:end))];

if rem(N_samples, 2) == 0
    full_scale = [amp_scale; flipud(amp_scale(2:end-1))];
else
    full_scale = [amp_scale; flipud(amp_scale(2:end))];
end

y = real(ifft(X .* full_scale));
y = y ./ std(y, 0, 1);
end