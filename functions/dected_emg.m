function [slopesChannelsxEpochs, other] = dected_emg(EEG)

% Select only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');

if ndims(EEG.data)==3
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
modulus = size(EEG.data(:,:),2)-N*L;
assert(modulus>=0);

% Compute (log) power spectra (7-75 Hz)
[NCHN,NPTS,NTRL]= size(dataeeg);
NWIN = NPTS;
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
for i = 1:NTRL
    [psdspectra(:,:,i),freq] = pwelch(dataeeg(:,:,i)',NWIN,0,NWIN,EEG.srate);
end

foi = [7 75];
frqmsk = freq>=foi(1) & freq<=foi(2);
psdspectra = permute(psdspectra,[1 3 2]);
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

end