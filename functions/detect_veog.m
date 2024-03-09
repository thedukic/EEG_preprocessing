function [dataeog, locs, treshold] = detect_veog(EEG)

chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataeog  = EEG.data(chaneog,:);

% Temporarily filter VEOG, bandpass 1-25 Hz
[bl, al] = butter(4,20/(EEG.srate/2),'low');
[bh, ah] = butter(4,1/(EEG.srate/2),'high');
assert(isstable(bl,al));
assert(isstable(bh,ah));
dataeog = filtfilt(bh,ah,filtfilt(bl,al,dataeog'))';

% % Temporarily filter VEOG, lowpass 6 Hz
% [bl, al] = butter(2,6/(EEG.srate/2),'low');
% assert(isstable(bl,al));
% dataeog = filtfilt(bl,al,dataeog);
% % dataeog = abs(dataeog);

% See: RELAX_blinks_IQR_method
EOGIQR = iqr(dataeog);
EOG75P = prctile(dataeog,75);
treshold = EOG75P + 4*EOGIQR;

% Approach 1: Detect eye blinks
% eyeBlinksEpochs = dataeog>=treshold;
%
% % Epoch into 1s
% L = EEG.srate;
% N = floor(length(eyeBlinksEpochs)/L);
% dataeog = reshape(dataeog(1:N*L),L,N);
% % eyeBlinksEpochs = any(reshape(eyeBlinksEpochs(1:N*L),L,N));
% eyeBlinksEpochs = reshape(eyeBlinksEpochs(1:N*L),L,N);
%
% % Binks should last at least this long [ms]
% mspersamp = 1000/EEG.srate;
% minpopdur = round(50/mspersamp);
% eyeBlinksEpochs = (find_maxduration(eyeBlinksEpochs')>minpopdur)';
%
% badtrl_all = find(eyeBlinksEpochs);
%
% fh = figure; hold on;
% F = (0:size(dataeog,1)-1)./EEG.srate;
% plot(F,dataeog(:,badtrl_all));

% Approach 2: Detect eye blinks
% times = EEG.times/1000; % EEGLAB time is in [ms] - > [s]
[qrspeaks, locs] = findpeaks(dataeog,EEG.times/1000,'MinPeakHeight',treshold,'MinPeakDistance',0.2);

end