function [psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG,thisScript)
% Helper function for estimating power spectra for:
% 1. preproc_cleaning2 script
% 2. report_final

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
chanemg = strcmp({EEG.chanlocs.type},'EMG');

if strcmpi(thisScript,'preproc2')
    if strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
        dataeeg = EEG.data;
    else
        dataeeg = EEG.data(chaneeg,:,:);
    end
elseif strcmpi(thisScript,'freport')
    % Make it continuous
    dataeeg = EEG.data(chaneeg,:);

    % Make it 4s long -> 0.25 Hz freq resolution
    [NCHN,NPTSALL] = size(dataeeg);
    NPTS = 4 * EEG.srate;
    NTRL = floor(NPTSALL/NPTS);
    dataeeg = reshape(dataeeg(:,1:NTRL*NPTS),NCHN,NPTS,NTRL);

end

% Compute power spectra
[NCHN, NPTS, NTRL] = size(dataeeg);
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

for i = 1:NTRL
    [psdspectra(:,:,i), freq] = pwelch(dataeeg(:,:,i)',NPTS,0,NPTS,EEG.srate);
end

% Average the spectra
psdspectra = mean(psdspectra,3);

end