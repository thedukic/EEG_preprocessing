function [eyeBlinksMask, eyeBlinksEpochs, BlinkMaxLatency, eyeBlinkData, treshold]= detect_veog(EEG,winBlink,ignoreExtremeNoise)
%
% SDukic, October 2024
%
% EOGfocus = +- time window around the max point of the blink [ms]
% EOGfocus = 400;  % ms for left/right eyeblink base
% EOGfocus = 2000; % ms for blink leftovers

% =========================================================================
% Apporach 1 (RELAX)
mspersmpl = 1000/EEG.srate;
winBlinksmpl = round(winBlink/mspersmpl);
fprintf('VEOG evaluation window is +-%d ms around the detected peaks.\n',winBlink);

chaneog = strcmp({EEG.chanlocs.labels},'VEOG');
eyeBlinkData = EEG.data(chaneog,:);

% Temporarily filter for better detection
% It is fine that it will be now double-filtered
[bl, al] = butter(4,15/(EEG.srate/2),'low');
[bh, ah] = butter(4,1/(EEG.srate/2),'high');

eyeBlinkData = do_filteringcore(bl,al,eyeBlinkData,EEG.event,EEG.srate);
eyeBlinkData = do_filteringcore(bh,ah,eyeBlinkData,EEG.event,EEG.srate);

% Is this a good idea tho?
% eyeBlinkData = abs(eyeBlinkData);

% Do this only if doing the MWF step
if ignoreExtremeNoise
    assert(length(eyeBlinkData)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1));
    eyeBlinkData(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = 0;
end

% VEOG
% figure; histogram(eyeBlinkData);
EOGIQR = iqr(eyeBlinkData);
EOG75P = prctile(eyeBlinkData,75);
if EOGIQR>100
    % If eyeblinks are too frequent and small?
    treshold = EOG75P + EOGIQR;
    fprintf('Eye blinks were detected using a treshold of 75PRC+IQR = %1.0fuV.\n',treshold);
else
    treshold = EOG75P + 2*EOGIQR;
    fprintf('Eye blinks were detected using a treshold of 75PRC+2IQR = %1.0fuV.\n',treshold);
end

% Treshold the EOG signal
BlinkIndexMetric = double(eyeBlinkData>treshold);

% figure; hold on;
% plot(EEG.times(1:40*256),dataeog(1:40*256));
% plot(EEG.times(1:40*256),treshold*BlinkIndexMetric(1:40*256));

% Check that blinks exceed the IQR threshold for more
% than 50ms, and to detect the blink peak within the period that
% exceeds the threshold:
ix_blinkstart = find(diff(BlinkIndexMetric)==1)+1;  % indices where BlinkIndexMetric goes from 0 to 1
ix_blinkend   = find(diff(BlinkIndexMetric)==-1);   % indices where BlinkIndexMetric goes from 1 to 0

if ~isempty(ix_blinkstart)
    if ix_blinkend(1,1)<ix_blinkstart(1,1); ix_blinkend(:,1)=[]; end % if the first downshift occurs before the upshift, remove the first value in end
    if ix_blinkend(1,size(ix_blinkend,2))<ix_blinkstart(1,size(ix_blinkstart,2)); ix_blinkstart(:,size(ix_blinkstart,2))=[];end % if the last upshift occurs after the last downshift, remove the last value in start

    BlinkThresholdExceededLength=ix_blinkend-ix_blinkstart; % length of consecutive samples where blink threshold was exceeded
    BlinkRunIndex = find(BlinkThresholdExceededLength>round(50/mspersmpl)); % find locations where blink threshold was exceeded by more than 50ms
    % find latency of the max voltage within each period where the blink
    % threshold was exceeded:
    if size(BlinkRunIndex,2)>0
        % continuousEEG.RELAX.IQRmethodDetectedBlinks = 1;
        % epochedEEG.RELAX.IQRmethodDetectedBlinks    = 1;
        for x=1:size(BlinkRunIndex,2)
            o=ix_blinkstart(BlinkRunIndex(x));
            c=ix_blinkend(BlinkRunIndex(x));
            [~,I] = max(eyeBlinkData(1,o:c),[],2);
            BlinkMaxLatency(1,x) = o+I;
        end
    end
end

eyeBlinksEpochs = [BlinkMaxLatency-winBlinksmpl; BlinkMaxLatency+winBlinksmpl]';

% Blinks cannot be outside of the recorded data
% Remove those that do not
eyeBlinksEpochs(eyeBlinksEpochs<1) = 1;
eyeBlinksEpochs(eyeBlinksEpochs>EEG.pnts) = EEG.pnts;

% Check if all blinks fit in the data
L = 2*winBlinksmpl+1;
NTRL = length(BlinkMaxLatency);
goodEyeBlinks = false(NTRL,1);

for i = 1:NTRL
    if length(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2))==L
        goodEyeBlinks(i) = true;
    end
end

eyeBlinksEpochs = eyeBlinksEpochs(goodEyeBlinks,:);

% Report
if ~all(goodEyeBlinks)
    warning('There are %d blinks excluded because their window was outside of the recorded data.',sum(~goodEyeBlinks));
end

% mask = zeros(size(dataeog));
% for i = 1:size(eyeBlinksEpochs,1)
%     mask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) = mask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2))+1;
% end
% figure; hold on;
% plot(EEG.times(274803:end),dataeog(274803:end));
% plot(EEG.times(274803:end),treshold*mask(274803:end));

% =========================================================================
% Approach 2 (old)
% chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
% dataeog  = EEG.data(chaneog,:);
% assert(sum(chaneog)==1);
%
% % Temporarily filter VEOG, bandpass 1-25 Hz
% [bl, al] = butter(2,25/(EEG.srate/2),'low');
% [bh, ah] = butter(2,0.5/(EEG.srate/2),'high');
% assert(isstable(bl,al));
% assert(isstable(bh,ah));
% dataeog = filtfilt(bh,ah,filtfilt(bl,al,dataeog'))';
%
% % See: RELAX_blinks_IQR_method
% EOGIQR = iqr(dataeog);
% EOG75P = prctile(dataeog,75);
% treshold = EOG75P + 2*EOGIQR;
%
% % mask = dataeog>treshold;
%
% % Detect eye blinks
% [qrspeaks, badLocs] = findpeaks(dataeog,EEG.srate,'MinPeakHeight',treshold,'MinPeakDistance',0.2);
%
% EOGfocus = 400; % ms
% EOGfocussamples = round(EOGfocus/(1000/EEG.srate));
% % fprintf('EOG duration is +-%dms wrt the detected peaks.\n',EOGfocus);
%
% eyeBlinksEpochs = round(badLocs*EEG.srate);
% eyeBlinksEpochs = [eyeBlinksEpochs-EOGfocussamples; eyeBlinksEpochs+EOGfocussamples]';
% eyeBlinksEpochs(eyeBlinksEpochs<1) = 1;
% eyeBlinksEpochs(eyeBlinksEpochs>length(dataeog)) = length(dataeog);

% =========================================================================

eyeBlinksMask = false(size(eyeBlinkData));
for i = 1:size(eyeBlinksEpochs,1)
    eyeBlinksMask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) = true;
end

% figure; hold on;
% % times = EEG.times/1000; % EEGLAB time is in [ms] - > [s]
% % plot(times(1:20*256),dataeog(1:20*256));
% timeaxis = (0:length(dataeog))/EEG.srate;
% plot(timeaxis(1:20*256),dataeog(1:20*256));
% plot(locs(1:11),qrspeaks(1:11),'ro');
%
% figure; hold on;
% plot(dataeog(1:30*256));
% plot(treshold*mask(1:30*256));

end