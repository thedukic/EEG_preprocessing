function [slopesChannelsxEpochs, other] = detect_emg(EEG,cfg)

% Select only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');

if ndims(EEG.data) == 3
    dataeeg = EEG.data(chaneeg,:,:);

    L = EEG.pnts;
    N = EEG.trials;
else
    dataeeg = EEG.data(chaneeg,:);

    % Epoch into 1s (or maybe better into 1s with 0.5 overlap, but OK)
    L = EEG.srate;
    N = floor(size(dataeeg,2)/L);
    dataeeg = reshape(dataeeg(:,1:N*L),sum(chaneeg),L,N);
end

% Check
modulus = size(EEG.data(:,:),2) - N*L;
assert(modulus >= 0);

% Compute (log) power spectra
[NCHN,NPTS,NTRL] = size(dataeeg);
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
for i = 1:NTRL
    [psdspectra(:,:,i),freq] = pwelch(dataeeg(:,:,i)',NPTS,0,NPTS,EEG.srate);
end

% The original reference uses 7-75 Hz and the RELAX toolbox as well, slope treshold >-0.31 or -0.51
% But MNE toolbox suggests 7-45 Hz and says that this is a better value in practice, but slope treshold is >0
% https://mne.tools/dev/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.find_bads_muscle
fprintf('Slope frequency range: %d-%d Hz\n',cfg.muscleSlopeFreq);
psdspectra = permute(psdspectra,[1 3 2]);

frqmsk = freq>=cfg.muscleSlopeFreq(1) & freq<=cfg.muscleSlopeFreq(2);
logpow = log10(psdspectra(frqmsk,:,:));
logfoi = log10(freq(frqmsk));

% Fit linear regression to log-log data, and store the slope
slopesChannelsxEpochs = NaN(NCHN,NTRL);
for i = 1:NCHN
    for j = 1:NTRL
        p = polyfit(logfoi,logpow(:,j,i),1);
        slopesChannelsxEpochs(i,j) = p(1);
    end
end

other.modulus = modulus;
other.L = L;
other.N = N;

% Exclude electrodes that are unlikely contamindated by EMG (?)
fprintf('Using only the slopes of peripheral electrodes.\n');
peripheralelecs = select_peripheralelecs(EEG);
slopesChannelsxEpochs(~peripheralelecs,:) = NaN;

end