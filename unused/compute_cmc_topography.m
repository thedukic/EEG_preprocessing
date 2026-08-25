function [cmc_topo, f_axis, cmc_all, cl95, fh] = compute_cmc_topography(DATA, cfg)
% COMPUTE_CMC_TOPOGRAPHY Calculates unrectified Corticomuscular Coherence (CMC)
% across all 128 scalp channels in the 15-45 Hz band using DPSS multitapers.
%
% Inputs:
%   DATA : EEGLAB dataset structure (continuous or epoched, 2048 Hz)
%   cfg  : (Optional) Configuration structure:
%          cfg.emg_chan     : 'APB' (E1-E2) or 'FDI' (E5-E6) [Default: 'APB']
%          cfg.event_marker : Cue trigger marker [Default: '21']
%          cfg.task_win_sec : Steady-state hold window in seconds [Default: [1.0, 5.0]]
%          cfg.freq_band    : Frequency range of interest [Default: [15, 45]]
%          cfg.emg_hp       : EMG high-pass filter in Hz [Default: 5]
%          cfg.apply_car    : Common Average Referencing across scalp leads [Default: true]
%          cfg.nw           : Multitaper time-halfbandwidth product [Default: 3.0]
%          cfg.num_tapers   : Number of Slepian tapers [Default: 5]
%          cfg.do_plot      : Render scalp topoplot [Default: true]
%
% Outputs:
%   cmc_topo : [128 x 1] vector of peak unrectified CMC (15-45 Hz) per channel
%   f_axis   : Frequency vector (Hz)
%   cmc_all  : [128 x n_freqs] full coherence matrix across all channels
%   cl95     : Theoretical 95% confidence limit
%   fh       : Figure handle

% -------------------------------------------------------------------------
% 1. Defaults and Setup
% -------------------------------------------------------------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'emg_chan'),     cfg.emg_chan     = 'FDI';        end
if ~isfield(cfg, 'event_marker'), cfg.event_marker = '21';         end
if ~isfield(cfg, 'task_win_sec'), cfg.task_win_sec = [2.0, 5.0];   end
if ~isfield(cfg, 'freq_band'),    cfg.freq_band    = [15, 45];     end
if ~isfield(cfg, 'emg_hp'),       cfg.emg_hp       = 5;            end
if ~isfield(cfg, 'apply_car'),    cfg.apply_car    = true;         end
if ~isfield(cfg, 'nw'),           cfg.nw           = 5.0;          end
if ~isfield(cfg, 'num_tapers'),   cfg.num_tapers   = floor(2 * cfg.nw - 1); end
if ~isfield(cfg, 'do_plot'),      cfg.do_plot      = true;         end

fs = DATA.srate;
all_labels = {DATA.chanlocs.labels};

% -------------------------------------------------------------------------
% 2. Coordinate Transfer for 128 Scalp Channels
% -------------------------------------------------------------------------
mat_data = load('biosemi128_eeglab.mat', 'chanlocs');
template_locs = mat_data.chanlocs;
template_labels = {template_locs.labels};
data_labels = {DATA.chanlocs.labels};

coord_fields = intersect(fieldnames(template_locs), ...
    {'X', 'Y', 'Z', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', ...
    'sph_theta_besa', 'sph_phi_besa', 'type', 'ref'});

scalp_ch_idx = zeros(1, 128);
for i = 1:128
    idx_match = find(strcmpi(data_labels, template_labels{i}), 1);
    if ~isempty(idx_match)
        scalp_ch_idx(i) = idx_match;
        for f = 1:length(coord_fields)
            fld = coord_fields{f};
            DATA.chanlocs(idx_match).(fld) = template_locs(i).(fld);
        end
        DATA.chanlocs(idx_match).type = 'EEG';
    end
end
assert(all(scalp_ch_idx > 0), 'Could not resolve all 128 BioSemi scalp channels.');

% -------------------------------------------------------------------------
% 3. Common Average Referencing (CAR) across Scalp Leads
% -------------------------------------------------------------------------
if cfg.apply_car
    car_mean = mean(DATA.data(scalp_ch_idx, :), 1);
    DATA.data(scalp_ch_idx, :) = DATA.data(scalp_ch_idx, :) - car_mean;
end

% -------------------------------------------------------------------------
% 4. Resolve Bipolar EMG Lead
% -------------------------------------------------------------------------
switch upper(cfg.emg_chan)
    case 'APB',          pairs = {'E1', 'E2'};  muscle_tag = 'APB';
    case {'FDI', 'FDRI'}, pairs = {'E5', 'E6'};  muscle_tag = 'FDI';
    otherwise
        error('compute_cmc:UnknownMuscle', 'Choose ''APB'' or ''FDI''.');
end

idx_pos = find(strcmpi(all_labels, pairs{1}) | strcmpi(all_labels, ['EXG' pairs{1}(2:end)]), 1);
idx_neg = find(strcmpi(all_labels, pairs{2}) | strcmpi(all_labels, ['EXG' pairs{2}(2:end)]), 1);
assert(~isempty(idx_pos) && ~isempty(idx_neg), 'Could not resolve EMG pair [%s, %s].', pairs{1}, pairs{2});

emg_raw_bipolar = DATA.data(idx_pos, :, :) - DATA.data(idx_neg, :, :);
emg_label = sprintf('%s (%s-%s)', muscle_tag, all_labels{idx_pos}, all_labels{idx_neg});

% -------------------------------------------------------------------------
% 5. Extract Steady-State Epochs (+1.0s to +5.0s)
% -------------------------------------------------------------------------
win_start_samp = round(cfg.task_win_sec(1) * fs);
win_end_samp   = round(cfg.task_win_sec(2) * fs);
epoch_pnts     = win_end_samp - win_start_samp + 1;

if ndims(DATA.data) == 3 && isfield(DATA, 'times')
    time_mask = (DATA.times >= cfg.task_win_sec(1) * 1000) & (DATA.times <= cfg.task_win_sec(2) * 1000);
    eeg_epochs = double(DATA.data(scalp_ch_idx, time_mask, :));
    emg_epochs = double(emg_raw_bipolar(1, time_mask, :));
else
    ev_types = {DATA.event.type};
    target_str = string(cfg.event_marker);
    ev_matches = false(1, length(ev_types));
    for ev_i = 1:length(ev_types)
        curr_ev = string(ev_types{ev_i});
        if curr_ev == target_str || strcmpi(curr_ev, "S " + target_str) || strcmpi(curr_ev, "S" + target_str)
            ev_matches(ev_i) = true;
        end
    end
    event_latencies = round([DATA.event(ev_matches).latency]);
    assert(~isempty(event_latencies), 'No events matching marker [%s] found.', target_str);

    n_total_pnts = size(DATA.data, 2);
    eeg_epochs = [];
    emg_epochs = [];

    for k = 1:length(event_latencies)
        t0 = event_latencies(k);
        idx_epoch = (t0 + win_start_samp) : (t0 + win_end_samp);
        if idx_epoch(1) >= 1 && idx_epoch(end) <= n_total_pnts
            eeg_epochs(:, :, end+1) = double(DATA.data(scalp_ch_idx, idx_epoch)); %#ok<AGROW>
            emg_epochs(:, :, end+1) = double(emg_raw_bipolar(1, idx_epoch));      %#ok<AGROW>
        end
    end
end

n_valid_trials = size(eeg_epochs, 3);
assert(n_valid_trials >= 3, 'Insufficient valid trials (N=%d).', n_valid_trials);

% -------------------------------------------------------------------------
% 6. EMG Conditioning (Unrectified)
% -------------------------------------------------------------------------
[b_hp, a_hp] = butter(4, cfg.emg_hp / (fs / 2), 'high');

% High-pass and demean unrectified EMG across samples
emg_filt = zeros(size(emg_epochs));
for tr = 1:n_valid_trials
    emg_filt(1, :, tr) = filtfilt(b_hp, a_hp, squeeze(emg_epochs(1, :, tr)));
end
emg_unrect = emg_filt - mean(emg_filt, 2);

% Demean EEG epochs across samples
eeg_sig = eeg_epochs - mean(eeg_epochs, 2);

% -------------------------------------------------------------------------
% 7. Vectorized DPSS Multitaper Estimation Across All 128 Channels
% -------------------------------------------------------------------------
N_pts = epoch_pnts;
[tapers, ~] = dpss(N_pts, cfg.nw, cfg.num_tapers);
K = cfg.num_tapers;

n_fft = max(2048, 2^nextpow2(N_pts));
n_freq_bins = (n_fft / 2) + 1;
f_axis = (0:(n_fft / 2))' * (fs / n_fft);

S_ee = zeros(128, n_freq_bins);
S_mm = zeros(1, n_freq_bins);
S_em = zeros(128, n_freq_bins);

for tr = 1:n_valid_trials
    eeg_tr = eeg_sig(:, :, tr);        % [128 x N_pts]
    emg_tr = emg_unrect(1, :, tr);      % [1 x N_pts]

    for k = 1:K
        w = tapers(:, k)';             % [1 x N_pts]

        % Vectorized FFT across all 128 channels simultaneously
        X = fft(eeg_tr .* w, n_fft, 2);
        Y = fft(emg_tr .* w, n_fft, 2);

        X = X(:, 1:n_freq_bins);       % [128 x n_freq_bins]
        Y = Y(1, 1:n_freq_bins);       % [1 x n_freq_bins]

        S_ee = S_ee + abs(X).^2;
        S_mm = S_mm + abs(Y).^2;
        S_em = S_em + (X .* conj(Y));
    end
end

% Compute Magnitude-Squared Coherence for all 128 channels
cmc_all = (abs(S_em).^2) ./ (S_ee .* repmat(S_mm, 128, 1)); % [128 x n_freq_bins]

% 95% Confidence Limit
total_L = n_valid_trials * K;
cl95 = 1 - (0.05)^(1 / (total_L - 1));

% -------------------------------------------------------------------------
% 8. Extract 15-45 Hz Peak Topography
% -------------------------------------------------------------------------
band_mask = (f_axis >= cfg.freq_band(1)) & (f_axis <= cfg.freq_band(2));
[cmc_topo, peak_freq_idx] = max(cmc_all(:, band_mask), [], 2);

f_band = f_axis(band_mask);
peak_freqs = f_band(peak_freq_idx);

[max_cmc_val, best_ch_idx] = max(cmc_topo);
best_ch_label = template_labels{best_ch_idx};
best_ch_freq  = peak_freqs(best_ch_idx);

fprintf('\n===================================================================\n');
fprintf('  Unrectified CMC Topography (%d-%d Hz | DPSS Multitaper)\n', cfg.freq_band(1), cfg.freq_band(2));
fprintf('===================================================================\n');
fprintf('  EMG Derivation:    %s (Unrectified)\n', emg_label);
fprintf('  Trials Analyzed:   N = %d (Total Estimates L = %d)\n', n_valid_trials, total_L);
fprintf('  95%% Conf. Limit:   CL95 = %.4f\n', cl95);
fprintf('  Top Hotspot:       Channel %s | Peak CMC = %.4f at %.1f Hz\n', ...
    best_ch_label, max_cmc_val, best_ch_freq);
fprintf('===================================================================\n\n');

% -------------------------------------------------------------------------
% 9. Visualisation
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot, return; end

fh = figure('Color', 'w', 'Position', [200, 200, 1000, 450], 'Name', 'Unrectified CMC Topography');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: Scalp Topoplot
ax1 = nexttile;
plot_title = sprintf('Unrectified CMC (%d-%d Hz)\nPeak: %s (%.3f @ %.1f Hz)', ...
    cfg.freq_band(1), cfg.freq_band(2), best_ch_label, max_cmc_val, best_ch_freq);
mytopoplot(cmc_topo(:), [], plot_title, ax1);

% Panel 2: Coherence Spectrum of the Peak Channel
nexttile;
plot_freq_mask = (f_axis >= 1) & (f_axis <= 50);
f_p = f_axis(plot_freq_mask);
cmc_best_p = cmc_all(best_ch_idx, plot_freq_mask);

y_lim_max = max([max(cmc_best_p), cl95 * 1.5, 0.05]) * 1.15;

patch([cfg.freq_band(1) cfg.freq_band(2) cfg.freq_band(2) cfg.freq_band(1)], ...
    [0 0 y_lim_max y_lim_max], [0.92 0.95 0.98], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(f_p, cmc_best_p, 'Color', [0.15 0.45 0.85], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Peak Lead: %s', best_ch_label));
yline(cl95, '--k', sprintf('CL_{95} = %.4f', cl95), 'LineWidth', 1.3);

xlim([1 50]);
ylim([0 y_lim_max]);
xlabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Magnitude-Squared Coherence', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Coherence Spectrum at Peak Channel [%s]', best_ch_label), 'FontSize', 11, 'FontWeight', 'bold');
grid on;
legend('Location', 'northeast', 'FontSize', 9);

end