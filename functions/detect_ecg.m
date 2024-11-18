function [ECGmask, ECGbadEpoch, ECGlatency, ECGdata, plusEstimate] = detect_ecg(EXT,winECG,ECGsignalLabel)
%
% Detect the QRS peaks in the ECG signal
% SDukic, October 2024
%

% Select ECG
chanecg = strcmp({EXT.chanlocs.labels},'ECG');

% Was ECG recorded?
if any(chanecg)
    ECGdata = EXT.data(chanecg,:);

    % MATLAB:
    % https://nl.mathworks.com/help/wavelet/ug/r-wave-detection-in-the-ecg.html
    wt = modwt(ECGdata,5);
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);
    y = imodwt(wtrec,'sym4');

    y = abs(y).^2;
    y = y./max(y);
    treshold = prctile(y,75) + 3*iqr(y);
    times = EXT.times/1000; % EEGLAB time is in [ms]
    [qrspeaks,ECGlatency] = findpeaks(y,times,'MinPeakHeight',treshold,'MinPeakDistance',0.8);

    % figure; hold on;
    % T = 10;
    % timeTmp = 0:1/256:T-1/256;
    % plot(timeTmp,y(1:T*256));
    % plot(timeTmp,treshold*ones(size(timeTmp)));

    if length(winECG) == 1
        % fprintf('ECG duration is +-%dms wrt the detected peaks.\n',winECG);
        winECG = winECG*[-1 1];
    end
    assert(length(winECG)==2);

    winECGsamples = round(abs(winECG)./(1000/EXT.srate));
    fprintf('ECG duration is %d-%d wrt the detected peaks.\n',winECG);

    ECGbadEpoch = round(ECGlatency*EXT.srate);
    ECGbadEpoch = [ECGbadEpoch-winECGsamples(1); ECGbadEpoch+winECGsamples(2)]';
    N = length(ECGbadEpoch(1,1):ECGbadEpoch(1,2));

    ECGbadEpoch(ECGbadEpoch<1) = 1;
    ECGbadEpoch(ECGbadEpoch>EXT.pnts) = EXT.pnts;

    ECG = NaN(size(ECGbadEpoch,1),N);
    theseExcl = [];
    for i = 1:size(ECGbadEpoch,1)
        if length(ECGbadEpoch(i,1):ECGbadEpoch(i,2))==N
            ECG(i,:) = ECGdata(ECGbadEpoch(i,1):ECGbadEpoch(i,2));
        else
            theseExcl = [theseExcl i];
        end
    end

    ECGbadEpoch(theseExcl,:) = [];
    ECG(theseExcl,:) = [];
    ECGlatency(theseExcl) = [];

    % Average ECG
    ECG  = ECG - mean(ECG,2);
    mECG = mean(ECG,1);

    % Detect noise
    cdist = 1-pdist2(ECG,mECG,'correlation');

    corrTrashold = prctile(cdist,85);
    if strcmpi(ECGsignalLabel,'Estimated')
        corrTrashold = max(min(corrTrashold, 0.98), 0.6);
    elseif strcmpi(ECGsignalLabel,'Recorded')
        corrTrashold = max(min(corrTrashold, 0.98), 0.85);
    else
        error('Something went wrong.');
    end

    theseExcl = cdist<corrTrashold;
    ECG(theseExcl,:) = [];
    ECGlatency(theseExcl) = [];
    fprintf('ECG correlation treshold %1.2f removed %d detections.\n',corrTrashold, sum(theseExcl));

    % This should not be the case if real ECG signal was used
    % (not estimated from EEG data)
    NECG = size(ECG,1);
    if NECG==0
        warning('Low quality ECG data... Quitting...');
        ECGmask     = NaN;
        ECGbadEpoch = NaN;
        ECGlatency  = NaN;
        ECGdata     = NaN;
        plusEstimate = NaN;
        return;
    end

    % figure; scatter(ECGlatency,ones(size(ECGlatency)));
    % figure; histogram(diff(ECGlatency));
    plusEstimate = diff(ECGlatency);
    plusEstimate(plusEstimate<0.5 | plusEstimate>1.5) = [];
    plusEstimate = mean(plusEstimate)*60;
    fprintf('Th final number of ECG events detected: %d (average pulse %1.1f per min)\n',NECG,plusEstimate);

    % Average ECG
    mECG = mean(ECG,1);

    fh = figure; hold on;
    F = (0:size(ECG,2)-1)./EXT.srate;
    F = F-mean(F);
    plot(F,ECG,'LineWidth',1.2);
    set(gca, 'ColorOrder',brewermap(NECG,'BuGn'));
    plot(F,mECG,'Color',[0.8 0.1 0.1],'LineWidth',3);
    title([ECGsignalLabel ' ECG signal, N = ' num2str(NECG)]);
    xlabel('Time (s)'); ylabel('ECG amplitude');

    % Save
    plotX=35; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EXT.ALSUTRECHT.subject.preproc,[EXT.ALSUTRECHT.subject.id '_detectedECGQRS']),'-dtiff','-r400');
    close(fh);

    ECGmask = false(size(ECGdata));
    for i = 1:size(ECGbadEpoch,1)
        ECGmask(ECGbadEpoch(i,1):ECGbadEpoch(i,2)) = true;
    end

else
    ECGmask      = NaN;
    ECGbadEpoch  = NaN;
    ECGlatency   = NaN;
    ECGdata      = NaN;
    plusEstimate = NaN;
end