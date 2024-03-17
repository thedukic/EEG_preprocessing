function [data, events] = detect_eog(EEG)
%
% Using https://github.com/bwrc/eogert
% SDukic, March 2024

% =========================================================================
% Apporach 2 (new)

% Select only EEG + EOG
chanheog = strcmp({EEG.chanlocs.labels},'HEOG');
chanveog = strcmp({EEG.chanlocs.labels},'VEOG');

data.heog = EEG.data(chanheog,:);
data.veog = EEG.data(chanveog,:);

% % Temporarily filter VEOG, bandpass 1-25 Hz
% [bl, al] = butter(2,25/(EEG.srate/2),'low');
% [bh, ah] = butter(2,0.1/(EEG.srate/2),'high');
% assert(isstable(bl,al));
% assert(isstable(bh,ah));
% tmp = filtfilt(bh,ah,filtfilt(bl,al,data.veog'))';
% 
% figure; hold on;
% times = EEG.times/1000;
% plot(times(:),tmp(:));
% figure; histogram(tmp);

% Detect VEOG/HEOG events
% This may fale if (not enough) blinks are detected in the training block
[BLI, SAC] = eogert_offline(data.heog,data.veog,EEG.srate,90);

% Did the detection fail?
if ~isempty(BLI)
    padding = 0.5; % in [s]

    % Blinks
    NVEOG = length(BLI.BLI_START);
    events.veog = cell(1,NVEOG);
    for i = 1:NVEOG
        events.veog{i} = [round((BLI.BLI_START(i)-padding)*EEG.srate), round((BLI.BLI_START(i)+BLI.BLI_DUR(i)+padding)*EEG.srate)];
    end

    % Saccades
    NHEOG = length(SAC.SAC_START);
    events.heog = cell(1,NHEOG);
    for i = 1:NHEOG
        events.heog{i} = [round((SAC.SAC_START(i)-padding)*EEG.srate), round((SAC.SAC_START(i)+SAC.SAC_DUR(i)+padding)*EEG.srate)];
    end
else
    events = [];
end

% % Make a mask
% noiseMask = false(1,EEG.pnts);
% for i = 1:NHEOG
%     noiseMask(events.heog{i}(1):events.heog{i}(2)) = true;
% end
%
% chanveog  = strcmp({EEG.chanlocs.labels},'VEOG');
% chanheog  = strcmp({EEG.chanlocs.labels},'HEOG');
% dataveog  = EEG.data(chanveog,:);
% dataheog  = EEG.data(chanheog,:);
%
% % Temporarily filter VEOG, bandpass 1-25 Hz
% [bl, al] = butter(2,25/(EEG.srate/2),'low');
% [bh, ah] = butter(2,1/(EEG.srate/2),'high');
% assert(isstable(bl,al)); assert(isstable(bh,ah));
% EEG.data(chanveog,:) = filtfilt(bh,ah,filtfilt(bl,al,dataveog'))';
% EEG.data(chanheog,:) = filtfilt(bh,ah,filtfilt(bl,al,dataheog'))';
%
% % Visual check
% EEG1 = EEG;
% EEG1.data(:,~logical(noiseMask)) = 0;
% vis_artifacts(EEG,EEG1,'ChannelSubset',[65:96, find(chanveog), find(chanheog)]);

% =========================================================================
% % Approach 1 (old)
% chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
% dataeog  = EEG.data(chaneog,:);
%
% % Temporarily filter VEOG, bandpass 1-25 Hz
% [bl, al] = butter(2,25/(EEG.srate/2),'low');
% [bh, ah] = butter(2,1/(EEG.srate/2),'high');
% assert(isstable(bl,al));
% assert(isstable(bh,ah));
% dataeog = filtfilt(bh,ah,filtfilt(bl,al,dataeog'))';
%
% % See: RELAX_blinks_IQR_method
% EOGIQR = iqr(dataeog);
% EOG75P = prctile(dataeog,75);
% treshold = EOG75P + 5*EOGIQR;
%
% % Detect eye blinks
% [qrspeaks, locs] = findpeaks(dataeog,EEG.times/1000,'MinPeakHeight',treshold,'MinPeakDistance',0.2);
%
% % figure; hold on;
% % times = EEG.times/1000; % EEGLAB time is in [ms] - > [s]
% % plot(times(1:20*256),dataeog(1:20*256));
% % % plot(times(1:20*256),dataeog(1:20*256));
% % % plot(locs(1:11),qrspeaks(1:11),'ro');

end