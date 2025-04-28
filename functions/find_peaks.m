function [locs, psdPromsMean] = find_peaks(psdspectra,freq)

% psdPeak = mean(psdspectraNorm,2);
% psdPeak  = psdPeak ./ sum(psdPeak,1);
maskFreq      = freq < 40;
freqTmp       = freq(maskFreq);
% psdspectraTmp = psdspectra ./ sum(psdspectra,1);
% psdspectraTmp = mean(psdspectra(maskFreq,:),2);
psdspectraTmp = mean(psdspectra,2);
psdspectraTmp = psdspectraTmp ./ max(psdspectraTmp);
psdspectraTmp = psdspectraTmp(maskFreq);

% figure; hold on;
% plot(freqTmp,psdspectraTmp);
% psdspectraTmp = movmean(psdspectraTmp,2);
% plot(freqTmp,psdspectraTmp);

% Initial step
[qrspeaks, locs, ~, proms] = findpeaks(psdspectraTmp,freqTmp);
% Guided detection
psdPromsMean   = mean(proms);
% psdPromsMean   = quantile(proms,0.2);
% psdPromsFinal  = max(psdPromsAll(i),0.0005);
% psdPromsFinal  = min(psdPromsFinal,0.1);
psdPromsFinal    = 0.025;
[qrspeaks, locs] = findpeaks(psdspectraTmp,freqTmp,'MinPeakProminence',psdPromsFinal);

end