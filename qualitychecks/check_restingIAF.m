function EEG = check_restingIAF(EEG)

% select data, set remaining params
eegchans = strcmp({EEG.chanlocs.type},'EEG');
dataTmp  = double(EEG.data(eegchans,:));

cmin   = 3;         % minimum number of channel estimates required for cross-channel averages (for tutorial data: min == 1, max == 6)
fRange = [1 48];    % spectral range (set to filter passband)
w  = [7 13];        % alpha peak search window (Hz)
Fw = 11;            % SGF frame width (11 corresponding to a frequency span of ~2.69 Hz @ ~.24Hz frequency resolution)
k  = 5;             % SGF polynomial order

% run `restingIAF`
[pSpec.sums, pSpec.chans, f] = restingIAF(dataTmp, sum(eegchans), cmin, fRange, EEG.srate, w, Fw, k);
pSpec.sums.freq = f;

% Report
fprintf('Individual alpha frequency: %1.3f Hz\n',pSpec.sums.paf);

% % Plot
% figure; hold on;
% plot(f,pSpec.sums.muSpec);
% [a,b] = min(abs(f - pSpec.sums.paf));
% plot(f(b), pSpec.sums.muSpec(b), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
% text(f(b) + 1, pSpec.sums.muSpec(b), sprintf('Alpha Peak: %.2f Hz', pSpec.sums.paf), 'FontSize', 12, 'Color', 'r');

% Output
EEG.ALSUTRECHT.pSpec = pSpec;

end