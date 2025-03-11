function [psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG,thisScript)
% Helper function for estimating power spectra for:
% 1. preproc_cleaning2 script
% 2. report_final
%
% N.B.
% If data is already epoched, especially if with overlapping segments
% Do not concatinate it and epoch it *differently*
% This cases spikes/ringing in the spectra
%

chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
chanemg = strcmp({EEG(1).chanlocs.type},'EMG');

if strcmpi(thisScript,'preproc2')
    % @preproc_cleaning2
    % Data is already epoched
    assert(ndims(EEG.data) == 3);

    fs = EEG.srate;
    if strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
        dataeeg = EEG.data;
    else
        dataeeg = EEG.data(chaneeg,:,:);
    end

elseif strcmpi(thisScript,'freport')
    % @preproc_cleaning1
    % @report_final

    fs = EEG.srate;
    % if strcmpi(EEG.ALSUTRECHT.subject.task,'MMN') || strcmpi(EEG.ALSUTRECHT.subject.task,'SART')  || strcmpi(EEG.ALSUTRECHT.subject.task,'RS')
    if ndims(EEG.data) == 3
        % @report_final
        % Data is already epoched
        NPTS = size(EEG.data,2);
    else
        % @preproc_cleaning1 (report leftovers)
        % Data not epoched
        NPTS = 2 * fs; % [samples]
    end

    % Make it continuous
    dataeeg = EEG.data - mean(EEG.data,2);
    dataeeg = dataeeg(chaneeg,:);

    [NCHN, NPTSALL] = size(dataeeg);
    NTRL = floor(NPTSALL/NPTS);
    dataeeg = reshape(dataeeg(:,1:NTRL*NPTS), NCHN,NPTS,NTRL);

elseif strcmpi(thisScript,'speaks')
    % preproc_cleaning1

    % Make it 20s long -> 0.05 Hz freq resolution
    fs = EEG(1).srate;
    winSizeCompleteSpectrum = 20; % [s]

    dataeeg = cat(2,EEG(:).data);
    dataeeg = dataeeg(chaneeg,:);
    NPTS = size(dataeeg,2);

    % We want at least 8 segments for proper usage of pwelch
    if winSizeCompleteSpectrum*fs > NPTS/8
        winSizeCompleteSpectrum = floor(NPTS/8/fs);
        warning('Dataset is short. Adjusted window size for whole data set spectrum calculation to be 1/8 of the length.')
    end

    % Make it 20s long -> 0.05 Hz freq resolution
    [NCHN, NPTSALL] = size(dataeeg);
    NPTS = winSizeCompleteSpectrum * fs;
    NTRL = floor(NPTSALL/NPTS);
    dataeeg = reshape(dataeeg(:,1:NTRL*NPTS),NCHN,NPTS,NTRL);

end

% Compute power spectra
[NCHN, NPTS, NTRL] = size(dataeeg);
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

dataeeg = double(dataeeg);
dataeeg = dataeeg - mean(dataeeg,2);
dataeeg = permute(dataeeg,[2 1 3]);

for i = 1:NTRL
    [psdspectra(:,:,i), freq] = pwelch(dataeeg(:,:,i), NPTS, 0, NPTS, fs);
end

% Average the spectra
psdspectra = mean(psdspectra,3);

end