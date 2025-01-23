

Ntmp = floor(size(ICAdata,2)/512);
dataTmp = reshape(ICAdata(:,1:512*Ntmp),70,512,Ntmp);

% Compute power spectra
[NCHN,NPTS,NTRL] = size(dataTmp);
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

for i = 1:NTRL
    [psdspectra(:,:,i), freq] = pwelch(dataTmp(:,:,i)',NPTS,0,NPTS,512);
end

psdspectra = mean(psdspectra,3);
figure; plot(freq,psdspectra(:,1));