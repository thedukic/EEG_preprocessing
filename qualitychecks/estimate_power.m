function [psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG,thisScript)
% Helper function for estimating power spectra for:
% 1. preproc_cleaning2 script
% 2. report_final

chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
chanemg = strcmp({EEG(1).chanlocs.type},'EMG');

if strcmpi(thisScript,'preproc2')
    if strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
        dataeeg = EEG.data;
    else
        dataeeg = EEG.data(chaneeg,:,:);
    end

elseif strcmpi(thisScript,'freport')
    % Make it 4s long -> 0.25 Hz freq resolution
    if strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
        % I think that because of the overlapping segments
        % If this is increased to longer trials, it causes spikes in spectra
        NPTS = size(EEG.data,2);
    else
        winSizeCompleteSpectrum = 4; % [s]
        NPTS = winSizeCompleteSpectrum * EEG.srate;
    end

    % Make it continuous
    dataeeg = EEG.data - mean(EEG.data,2);
    dataeeg = dataeeg(chaneeg,:);

    [NCHN, NPTSALL] = size(dataeeg);
    NTRL = floor(NPTSALL/NPTS);
    dataeeg = reshape(dataeeg(:,1:NTRL*NPTS), NCHN,NPTS,NTRL);

elseif strcmpi(thisScript,'speaks')
    % Make it 20s long -> 0.05 Hz freq resolution
    winSizeCompleteSpectrum = 20; % [s]

    dataeeg = cat(2,EEG(:).data);
    dataeeg = dataeeg(chaneeg,:);
    NPTS = size(dataeeg,2);

    % we want at least 8 segments for proper usage of pwelch
    if winSizeCompleteSpectrum*EEG(1).srate > NPTS/8
        winSizeCompleteSpectrum = floor(NPTS/8/EEG(1).srate);
        warning('Dataset is short. Adjusted window size for whole data set spectrum calculation to be 1/8 of the length.')
    end

    % Make it 20s long -> 0.05 Hz freq resolution
    [NCHN, NPTSALL] = size(dataeeg);
    NPTS = winSizeCompleteSpectrum * EEG(1).srate;
    NTRL = floor(NPTSALL/NPTS);
    dataeeg = reshape(dataeeg(:,1:NTRL*NPTS),NCHN,NPTS,NTRL);

end

% Compute power spectra
[NCHN, NPTS, NTRL] = size(dataeeg);
% psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
psdspectra = NaN(NPTS/2+1,NCHN,NTRL);

dataeeg = double(dataeeg);
dataeeg = dataeeg - mean(dataeeg,2);
dataeeg = permute(dataeeg,[2 1 3]);

for i = 1:NTRL
    [psdspectra(:,:,i), freq] = pwelch(dataeeg(:,:,i),NPTS,0,NPTS,EEG(1).srate);
end

% Average the spectra
psdspectra = mean(psdspectra,3);

end