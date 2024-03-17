function [EEG, noiseMask] = mwf_heartbeats(EEG,theseData)
%
% Test whether only eye blinks are present using bipolar VEOG
% SDukic, March 2024
%

fprintf('\nMWF (ECG) heart beats...\n');

% Select ECG
chanecg = strcmp({EEG.chanlocs.labels},'ECG');

% Was ECG recorded?
if any(chanecg)
    % Select only EEG/EMG
    if strcmpi(theseData,'EEG')
        chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
        loglabel = 'R2';
    elseif strcmpi(theseData,'EMG')
        chaneeg  = strcmp({EEG.chanlocs.type},'EMG');
        loglabel = 'EMG';
    end

    % ECG recorded
    dataeeg  = EEG.data(chaneeg,:);
    dataecg  = EEG.data(chanecg,:);

    % MATLAB:
    % https://nl.mathworks.com/help/wavelet/ug/r-wave-detection-in-the-ecg.html
    wt = modwt(dataecg,5);
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);
    y = imodwt(wtrec,'sym4');

    y = abs(y).^2;
    y = y./max(y);
    treshold = prctile(y,75) + 3*iqr(y);
    times = EEG.times/1000; % EEGLAB time is in [ms]
    [qrspeaks,locs] = findpeaks(y,times,'MinPeakHeight',treshold,'MinPeakDistance',0.8);

    ECGfocus = 80; % ms
    ECGfocussamples = round(ECGfocus/(1000/EEG.srate));
    fprintf('ECG duration is +-%dms wrt the detected peaks.\n',ECGfocus);

    badEpoch2 = locs*EEG.srate;
    badEpoch2 = [badEpoch2-ECGfocussamples; badEpoch2+ECGfocussamples]';
    N = length(badEpoch2(1,1):badEpoch2(1,2));

    badEpoch2(badEpoch2<1) = 1;
    badEpoch2(badEpoch2>EEG.pnts) = EEG.pnts;

    ECG = NaN(size(badEpoch2,1),N);
    for i = 1:size(badEpoch2,1)
        if length(badEpoch2(i,1):badEpoch2(i,2))==N
            ECG(i,:) = dataecg(badEpoch2(i,1):badEpoch2(i,2));
        end
    end
    ECG(isnan(ECG(:,1)),:) = [];
    mECG = median(ECG,1);

    cdist = 1-pdist2(ECG,mECG,'correlation');
    ctrsh = max(prctile(cdist,75),0.85);
    fprintf('ECG correlation treshold: %1.2f\n',ctrsh);
    ECG(cdist<ctrsh,:) = [];
    NECG = size(ECG,1);

    fh = figure; hold on;
    F = (0:size(ECG,2)-1)./EEG.srate;
    F = F-mean(F);
    plot(F,ECG','LineWidth',1.2);
    set(gca, 'ColorOrder',brewermap(NECG,'BuGn'));
    plot(F,mECG,'Color',[0.8 0.1 0.1],'LineWidth',3);
    title(['N = ' num2str(NECG)]);
    xlabel('Time (s)'); ylabel('ECG amplitude');

    % Save
    plotX=35; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ECG_' theseData]),'-dtiff','-r400');
    close(fh);

    % badEpoch2 = locs*EEG.srate;
    % badEpoch2 = [badEpoch2-10; badEpoch2+15]';
    % N = length(badEpoch2(1,1):badEpoch2(1,2));
    % badEpoch2(badEpoch2<1) = 1;
    % badEpoch2(badEpoch2>EEG.pnts) = EEG.pnts;
    % TECG = (N-1)/2/EEG.srate*1000;
    % fprintf('ECG duration is +-%2.0fms wrt the detected peaks.\n',TECG);

    % myCmap = brewermap(size(ECG,1),'Purples');
    % F = (0:size(ECG,2)-1)./EEG.srate;
    % amean = mean(ECG,1);
    % N     = size(ECG,1);      % Number of particiapnts
    % astd  = std(ECG)/sqrt(N); % Compute Standard Error Of The Mean (SEM)
    % CI95  = tinv(0.975,N-1);  % Calculate 95% Probability Intervals of t-Distribution
    % astd  = CI95*astd;        % Calculate 95% Confidence Intervals of all particiapnts at each x
    % plot(F,ECG);
    % plot(F,median(ECG),'Color',[0 0 0],'LineWidth',1.1);
    % fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],'k','FaceAlpha',0.5,'linestyle','none','HandleVisibility','off');
    % plot(F,amean,'Color',[0.8 0.2 0.2],'linewidth',1.2);

    %% Multi-channel Wiener filter
    % Noise mask:
    % NaN - ignore these samples
    % 0   - good data
    % 1   - bad data

    % Noise mask
    noiseMask = zeros(1,EEG.pnts);
    if ~isempty(badEpoch2)
        for i = 1:size(badEpoch2,1)
            noiseMask(badEpoch2(i,1):badEpoch2(i,2)) = 1;
        end
    end

    if strcmpi(theseData,'EEG')
        assert(length(noiseMask)==size(EEG.data,2));
        assert(length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)==length(noiseMask));
        noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;
    end

    % Log
    EEG.ALSUTRECHT.MWF.(loglabel).badElectrodes          = NaN;
    EEG.ALSUTRECHT.MWF.(loglabel).noiseMask              = noiseMask;
    EEG.ALSUTRECHT.MWF.(loglabel).proportionMarkedForMWF = mean(noiseMask,'omitnan');

    fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.(loglabel).proportionMarkedForMWF);
    fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF %s (ECG) heart beats\n',theseData);
    fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'Heart beat amplitude threshold:   %1.2f\n',treshold);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Heart beat correlation threshold: %1.2f\n',ctrsh);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Number of heart beats detected:   %d\n',NECG);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Heart beat duration for MWF: %1.2f\n',2*ECGfocus);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.(loglabel).proportionMarkedForMWF);

    if EEG.ALSUTRECHT.MWF.(loglabel).proportionMarkedForMWF>0.05
        [cleanEEG, d, W, SER, ARR] = mwf_process(dataeeg,noiseMask,8);

        % Check
        EEG0 = EEG;
        EEG0.data(chaneeg,:) = cleanEEG;
        vis_artifacts(EEG0,EEG);

        EEG.data(chaneeg,:) = cleanEEG;

        % EEG.ALSUTRECHT.MWF.(loglabel).estimatedArtifactInEachChannel = d;
        % EEG.ALSUTRECHT.MWF.(loglabel).matrixUsedToEstimateArtifacts  = W;
        EEG.ALSUTRECHT.MWF.(loglabel).status                 = 1;
        EEG.ALSUTRECHT.MWF.(loglabel).signalToErrorRatio     = SER;
        EEG.ALSUTRECHT.MWF.(loglabel).artifactToResidueRatio = ARR;

        fprintf(EEG.ALSUTRECHT.subject.fid,'Signal to error ratio:     %1.2f\n',SER);
        fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact to residue ratio: %1.2f\n',ARR);
    else
        EEG.ALSUTRECHT.MWF.(loglabel).status                 = 0;
        EEG.ALSUTRECHT.MWF.(loglabel).signalToErrorRatio     = NaN;
        EEG.ALSUTRECHT.MWF.(loglabel).artifactToResidueRatio = NaN;

        fprintf(EEG.ALSUTRECHT.subject.fid,'MWF will not be done. Too little data.\n');
        fprintf('MWF will not be done. Too little data.\n');
    end
else
    fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (ECG) heart beats\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF is not done as this participant does not have ECG recorded (but likely L/R earlobes instead).\n');
    fprintf('MWF is not done as this participant does not have ECG recorded (but likely L/R earlobes instead).\n');

    EEG.ALSUTRECHT.MWF.(loglabel).status                 = 0;
    EEG.ALSUTRECHT.MWF.(loglabel).badElectrodes          = NaN;
    EEG.ALSUTRECHT.MWF.(loglabel).noiseMask              = NaN;
    EEG.ALSUTRECHT.MWF.(loglabel).proportionMarkedForMWF = NaN;
    EEG.ALSUTRECHT.MWF.(loglabel).signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.(loglabel).artifactToResidueRatio = NaN;
end

end