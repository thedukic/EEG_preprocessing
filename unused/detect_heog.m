function [saccadesMask, saccadesEpochs, BlinkMaxLatency, saccadesData, treshold] = detect_heog(EEG,winSaccades)

BlinkMaxLatency = NaN;

% Demean and take absolute value so that L/R HEOG can be both detected (?)
chanheog = strcmp({EEG.chanlocs.labels},'HEOG');
assert(any(chanheog));

saccadesData = EEG.data(chanheog,:);
% dataheog = dataheog - trimmean(dataheog,10);
% dataheog = abs(dataheog);

% Temporarily filter HEOG, bandpass 1-25 Hz
[bl, al] = butter(2,20/(EEG.srate/2),'low');
[bh, ah] = butter(2,1/(EEG.srate/2),'high');

assert(isstable(bl,al));
assert(isstable(bh,ah));

saccadesData = do_filteringcore(bl,al,saccadesData,EEG.event,EEG.srate);
saccadesData = do_filteringcore(bh,ah,saccadesData,EEG.event,EEG.srate);

% dataheog = dataheog - trimmean(dataheog,10);

% HEOG treshold
EOGIQR = iqr(saccadesData);
EOG75P = prctile(saccadesData,75);
treshold = EOG75P + 1.5*EOGIQR;
saccadesEpochs = saccadesData>=treshold;

% Needed to exclude false positives that are actually eye blinks
[noiseMaskBlink, eyeBlinksEpochs] = detect_veog(EEG,200,false);

saccadesMask = zeros(1,EEG.pnts);
if any(saccadesEpochs)
    jump = find(diff([false, saccadesEpochs, false])~=0);
    durall = jump(2:2:end)-jump(1:2:end);

    % Minimum of 25 ms of HEOG
    mspersamp = 1000/EEG.srate;
    mindrftdur = round(25/mspersamp);

    % Horizontal eye movements should last at least this long [ms]
    saccadesEpochs = find(durall>mindrftdur);
    NHEOG = length(saccadesEpochs);

    if NHEOG>0
        jumpStart = jump(1:2:end);
        jumpStop  = jump(2:2:end);

        % winSaccades = 300; % ms
        EOGfocussamples = round(winSaccades/mspersamp);

        jumpStart = jumpStart-EOGfocussamples;
        jumpStart(jumpStart<1) = 1;
        jumpStop = jumpStop+EOGfocussamples;
        jumpStop(jumpStop>EEG.pnts) = EEG.pnts;

        actuallyBlinks = false(NHEOG,1);
        for i = 1:NHEOG
            actuallyBlinks(i) = any(noiseMaskBlink(jumpStart(saccadesEpochs(i)):jumpStop(saccadesEpochs(i))));
        end

        saccadesEpochs(actuallyBlinks) = [];
        % jumpStart(actuallyBlinks)       = [];
        % jumpStop(actuallyBlinks)        = [];
        NHEOG = length(saccadesEpochs);

        % % Check: A lot of times HEOG captures blinks too !!!
        % EEGTMP = EEG;
        % mask = false(size(EEGTMP.times));
        % for i = 1:NHEOG
        %     mask(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i))) = true;
        % end
        % EEGTMP.data(:,~mask) = 0;
        % vis_artifacts(EEG,EEGTMP);

        fh = figure; hold on;
        durheog = NaN(NHEOG,1);
        for i = 1:NHEOG
            y = saccadesData(jumpStart(saccadesEpochs(i)):jumpStop(saccadesEpochs(i)));
            durheog(i) = length(y);
            T = linspace(-0.5,0.5,durheog(i));

            absDiff_100ms = abs(T - (-0.4));
            minDiff_100ms = min(absDiff_100ms(:));
            [~, col_m150ms] = find(absDiff_100ms == minDiff_100ms);
            absDiff_100ms = abs(T - 0.4);
            minDiff_100ms = min(absDiff_100ms(:));
            [~, col_p100ms] = find(absDiff_100ms == minDiff_100ms);
            y = y - mean(y([1:col_m150ms, col_p100ms:end]));
            % y = y - mean(y(1:col_m150ms));

            plot(T,y,'LineWidth',1.2);
        end
        set(gca,'ColorOrder',brewermap(NHEOG,'BuGn'));
        title(['N = ' num2str(NHEOG)]);
        xlabel('Time a.u.'); ylabel('HEOG amplitude');

        % Save
        plotX=35; plotY=20;
        set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_detectedHEOG']),'-dtiff','-r400');
        close(fh);

        TEOG = mean(durheog)/EEG.srate*1000;
        fprintf('HEOG average duration is %2.0fms.\n',TEOG);

        for i = 1:NHEOG
            saccadesMask(jumpStart(saccadesEpochs(i)):jumpStop(saccadesEpochs(i))) = 1;
        end
    end

end