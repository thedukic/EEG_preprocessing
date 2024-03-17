function [EEG, badElectrodes, noiseMask] = detect_channelpops(EEG,cfgbch)
% Test whether electrodes have high voltage jumps or strong wobbles
% Done using simple tresholding
%
fprintf('\nDetecting and fixing large electrode pops...\n');

% Some channels are external (eg ECG)
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chanlocs = EEG.chanlocs(chaneeg);
labels   = {EEG.chanlocs.labels};
chaneog  = strcmp(labels,'VEOG');

% Select only EEG + VEOG
dataeog = EEG.data(chaneog,:);
dataeeg = EEG.data(chaneeg,:);
NCHNEEG = sum(chaneeg);

if strcmpi(cfgbch.popType,'large')
    % Filter EOG, bandpass 1-25 Hz
    FNYQ = EEG.srate/2;
    [bl, al] = butter(2,7/FNYQ,'low');
    assert(isstable(bl,al));
    dataeog = filtfilt(bl,al,dataeog);

    % Use frontal electrodes for blink detection ???
    % BlinkElectrodes = {'C15','C16','C17','C18','C28','C29'};
    % BlinkElectrodes = labels(contains(labels,BlinkElectrodes));
    % EEGEyeOnly.data = mean(EEGEyeOnly.data,1); % This could be changed to median, which would increase robustness against bad channels or outliers (but wouldn't work so well if including electrodes barely affected by blinks)

    % See: RELAX_blinks_IQR_method
    EOGIQR = iqr(dataeog,2);
    EOG75P = prctile(dataeog,75,2);
    UpperBound = EOG75P + 3*EOGIQR;
    BlinkIndexMetric = dataeog>=UpperBound;

elseif strcmpi(cfgbch.popType,'small')
    dataeog = 0*dataeog;
    BlinkIndexMetric = true(size(dataeog));
end

% Baseline correct (robust) EEG only
dataeeg = dataeeg - trimmean(dataeeg,10,'round',2);

% Epoch into 1s
L = EEG.srate;
N = floor(size(dataeeg,2)/L);
dataeeg = reshape(dataeeg(:,1:N*L),NCHNEEG,L,N);
dataeog = reshape(dataeog(1:N*L),L,N);
BlinkIndexMetric = any(reshape(BlinkIndexMetric(1:N*L),L,N));
time = reshape(EEG.times(1:N*L),L,N)/1000; % [ms] -> [s]

assert(size(dataeeg,3)==size(dataeog,2));
assert(size(dataeeg,3)==length(BlinkIndexMetric));

% Define and initialise
cnt    = 0;
NTRL   = size(dataeeg,3);
badchn = zeros(NCHNEEG,1);
badtrl_rmv = [];
badtrl_fix = [];

if strcmpi(cfgbch.popType,'large')
    popThreshold     = cfgbch.popThresholdLarge;
    popThresholdPlot = popThreshold;
    K = 0; % K = [-1 0 1];
elseif strcmpi(cfgbch.popType,'small')
    popThreshold     = cfgbch.popThresholdSmall;
    popThresholdPlot = 50;
    K = 0;
else
    error('Wrong settings for electrode pop detection!');
end

% Plot
fh = figure;
th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

% Drifts/pops should last at least this long [ms]
mspersamp = 1000/EEG.srate;
minpopdur = round(50/mspersamp);

% Which epochs are without any eye blinks
BlinkIndexMetric = ~BlinkIndexMetric;

for i = 1:NTRL
    if BlinkIndexMetric(i)
        mask = abs(dataeeg(:,:,i)) >= popThreshold;
        if any(mask,"all")
            if strcmpi(cfgbch.popType,'large')
                % Ichn = any(abs(dataeeg(:,:,i)) >= popThreshold,2);
                % Ichn = sum(abs(dataeeg(:,:,i)) >= popThreshold,2)>minpopdur;
                Ichn = find_maxduration(mask)>minpopdur;

            elseif strcmpi(cfgbch.popType,'small')
                % m0 = abs(robust_zscore(max(dataeeg(:,:,i),[],2) - min(dataeeg(:,:,i),[],2),[])) >= popThreshold;
                % m1 = abs(robust_zscore(dataeeg(:,:,i),[])) >= popThreshold;
                % m1 = any(m1,2);
                % m2 = abs(robust_zscore(diff(dataeeg(:,:,i)'),[]))' >= popThreshold;
                % m2 = any(m2,2);
                m1 = abs(robust_zscore(diff(dataeeg(:,:,i)')',[])) >= popThreshold;
                m1 = any(m1,2);
                % x = abs(corr(dataeeg(:,:,i)',type="Spearman"));
                % x = x-diag(diag(x));
                % m2 = quantile(x, 0.98)'<0.6;
                % Ichn = [m0|m1|m2; false];
                Ichn = m1;
            end

            if any(Ichn)
                % Not sure if the below line should be left here or within
                % the if-else below; if here then the code might be too strict
                % as it counts in all large deviations regardless if a (very bad) epoch
                % will be removed from further analysis ...
                badchn(Ichn) = badchn(Ichn)+1;

                if sum(Ichn)>2
                    % badchn(Ichn) = badchn(Ichn)+1;
                    badtrl_rmv = [badtrl_rmv, i+K];
                else
                    badtrl_fix = [badtrl_fix, i+K];
                    % badchn(Ichn) = badchn(Ichn)+1;
                    cnt = cnt+1;

                    % Red:   Bad channels
                    % Black: EOG channel
                    % Green: Corrected channel
                    nexttile; hold on;
                    plot(time(:,i),dataeeg(Ichn,:,i),'Color',[0.6 0.1 0.1],'LineWidth',2);
                    plot(time(:,i),dataeog(:,i),'Color',[0.1 0.1 0.1],'LineWidth',1.1);
                    % plot(tmp.time{i},tmp.trial{i}(~any(m,2),:),'Color',[0.1 0.1 0.1],'LineWidth',1.1);
                    axis tight; ylim(popThresholdPlot.*[-1 1]);
                    % title([EEGTMP.label{Ichn} ', ' num2str(mean(avgchange))]);
                    % title(num2str(length(EEGTMP.label(Ichn))));
                    title([labels{Ichn}]);

                    l = {find(Ichn)};
                    if length(K)>1
                        % Interpolate three trials: i-1, i, i+1 trials
                        for j = 1:length(K)
                            k = i+K(j);
                            if k>0 && k<=NTRL
                                TMP = interpolate_epochs(dataeeg(:,:,k),chanlocs,l);
                                % dataeeg(:,:,k) = [TMP; dataeeg(NCHNEOG,:,k)];
                                dataeeg(:,:,k) = TMP;
                            end
                        end
                    else
                        % Interpolate one trial: i trial
                        TMP = interpolate_epochs(dataeeg(:,:,i),chanlocs,l);
                        % dataeeg(:,:,i) = [TMP; dataeeg(NCHNEOG,:,i)];
                        dataeeg(:,:,i) = TMP;
                    end

                    plot(time(:,i),dataeeg(Ichn,:,i),'Color',[0.1 0.6 0.1],'LineWidth',1.1);
                    % xticks(mean(round(time(:,i)([1 end]))));
                    xticks(round(time([1 end],i)));
                    yticks(popThresholdPlot*[-1 1]); ylim(1.5*popThresholdPlot*[-1 1]);
                end
            end
        end
    end
end

% Save
plotX=35; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_detectedpops_' cfgbch.popType '_1']),'-dtiff','-r400');
close(fh);

% Trials can be overlapping
badtrl_fix = unique(badtrl_fix);
badtrl_rmv = unique(badtrl_rmv);

% Determine bad channels based on the amount of pops
NBAD = length(badtrl_rmv);
badchn = badchn./NTRL;
badElectrodes = {EEG.chanlocs(badchn>cfgbch.popTime).labels};
% badElectrodes = {EEG.chanlocs(badchn>=cfgbch.popTime).labels};

% Report
fprintf('Simple channel pops were detected: %d time(s)\n',cnt);
fprintf('Large channel deviations accross many channels were detected: %d time(s)\n',NBAD);
fprintf('Maximum time of detected pops: %1.3f\n',max(badchn));

% Log
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Electrode pops (%s)\n',cfgbch.popType);
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye blink threshold: %1.2f\n',UpperBound);
fprintf(EEG.ALSUTRECHT.subject.fid,'Percentage od epochs with eye blinks detected: %1.2f\n',1-mean(BlinkIndexMetric));
fprintf(EEG.ALSUTRECHT.subject.fid,'Simple channel pops were detected: %d time(s)\n',cnt);
fprintf(EEG.ALSUTRECHT.subject.fid,'Large channel deviations accross many channels were detected: %d time(s)\n',NBAD);
fprintf(EEG.ALSUTRECHT.subject.fid,'Maximum time of detected pops: %1.3f\n',max(badchn));

% Plot
if NBAD>0
    % Plot max 19 trials that were marked as very bad
    NBAD2 = min(NBAD,19);

    % Plot
    fh = figure;
    th = tiledlayout('flow');
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    for i = 1:NBAD2
        nexttile;
        plot(time(:,i),dataeeg(:,:,badtrl_rmv(i)),'Color',[0.5 0.5 0.5],'LineWidth',1.1);
        xticks(round(time([1 end],i)));
        yticks(popThresholdPlot*[-1 1]); ylim(1.5*popThresholdPlot*[-1 1]);
    end

    if NBAD>i
        nexttile;
        plot(time(:,i),randn(1,length(time(:,i))),'Color',[0.9 0.5 0.5],'LineWidth',1.1);
        xticks(round(time([1 end],i)));
        yticks(popThresholdPlot*[-1 1]); ylim(1.5*popThresholdPlot*[-1 1]);
    end

    % Save
    plotX=35; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_detectedpops_' cfgbch.popType '_2']),'-dtiff','-r400');
    close(fh);
end

% Put interpolated segments
if NTRL*EEG.srate==EEG.pnts
    dataeeg = dataeeg(:,:);
else
    n = NTRL*EEG.srate+1:EEG.pnts;
    dataeeg = [dataeeg(:,:), EEG.data(chaneeg,n)];
end
assert(size(dataeeg,1)==sum(chaneeg));
assert(size(dataeeg,2)==EEG.pnts);

% Put back to the original dataset, but not the filtered EOG
EEG.data(chaneeg,:) = dataeeg;
clearvars data;

% Double-checks
assert(size(EEG.data,1)==EEG.nbchan);
assert(size(EEG.data,2)==EEG.pnts);
EEG = eeg_checkset(EEG);

% Make a noise mask using very bad epochs
badtrl_all = badtrl_rmv;
noiseMask  = zeros(1,EEG.pnts);
if ~isempty(badtrl_all)
    badEpoch1 = [badtrl_all-1; badtrl_all]'; % in [s]
    badEpoch2 = badEpoch1*EEG.srate;         % in [samples]
    badEpoch2(:,1) = badEpoch2(:,1)+1;
    if ~isempty(badtrl_all)
        for i = 1:NBAD
            noiseMask(badEpoch2(i,1):badEpoch2(i,2)) = 1;
        end
    end
    assert(length(noiseMask)==size(EEG.data,2));
end

% =========================================================================
%  Save and update the extreme noise mask
if strcmpi(cfgbch.popType,'large')
    EEG.ALSUTRECHT.pops.large.badElectrodes = badElectrodes;
    EEG.ALSUTRECHT.pops.large.noiseMask     = noiseMask;
else
    EEG.ALSUTRECHT.pops.small.badElectrodes = badElectrodes;
    EEG.ALSUTRECHT.pops.small.noiseMask     = noiseMask;
end

if strcmpi(cfgbch.popType,'large')
    % Some epochs with extreme pops may be fixed here
    if ~isempty(badtrl_fix)
        [i,j] = find(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3==badtrl_fix');
        if ~isempty(j)
            EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3(j) = [];
            fprintf(EEG.ALSUTRECHT.subject.fid,'Some epochs with extreme pops are fixed: %d\n',length(j));
        end
    end

    % Add possibly newly detected very bad epochs
    if ~isempty(badtrl_rmv)
        NTMP = length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
        EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3 = unique([EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3, badtrl_rmv]);
        % EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs = unique([EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs; badepoch1],'rows');

        if NTMP<length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
            fprintf(EEG.ALSUTRECHT.subject.fid,'Some new epochs with extreme pops are added: %d\n',length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)-NTMP);
        end
    end

    % Update the extreme noise mask only if something was detected
    if ~isempty(badtrl_fix) || ~isempty(badtrl_rmv)
        extremeNoiseMask1 = false(1,EEG.pnts);
        extremeEpochs = EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3;
        if ~isempty(extremeEpochs)
            badEpoch1 = [extremeEpochs-1; extremeEpochs]'; % in [s]
            badEpoch2 = badEpoch1*EEG.srate;               % in [samples]
            badEpoch2(:,1) = badEpoch2(:,1)+1;

            for i = 1:length(extremeEpochs)
                extremeNoiseMask1(badEpoch2(i,1):badEpoch2(i,2)) = true;
            end
        end

        extremeNoiseMask2 = false(1,NTRL);
        extremeNoiseMask2(extremeEpochs) = true;

        % Update these too
        EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier = mean(extremeNoiseMask1);
        EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1                 = extremeNoiseMask1;
        EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2                 = extremeNoiseMask2;
        EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs4                 = badEpoch2;
    end

    % Double-check
    assert(sum(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)/EEG.srate==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3));
    assert(sum(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3));
    assert(size(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs4,1)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3));

end

%%
% =========================================================================
%                              HELPER FUNCTIONS
% =========================================================================
% 1
function y = robust_zscore(x1,x2)
% [] W.A. Stahel. Robuste Schatzungen: infinitesimale Optimalit¨at und Schatzungen von Kovarianzmatrizen. PhD thesis, ETH Zurich, 1981.
% [] D.L. Donoho. Breakdown properties of multivariate location estimators. Qualifying paper, Harvard University, Boston, 1982.

% Operates along columns if the inputs are matrices
if ~isvarname(x2)
    x2 = [];
end

if isempty(x2)
    % disp('Robust Z-score using patient data.');
    y = (x1-median(x1)) ./ (1.4826*mad(x1,1));     % median absolute deviation
    % y = (x1-median(x1)) ./(1.253314*mad(x1,0));  % mean absolute deviation
else
    % disp('Robust Z-score using control data.');
    y = (x1-median(x2)) ./ (1.4826*mad(x2,1));      % median absolute deviation
    % y = (x1-median(x2)) ./ (1.253314*mad(x2,0));  % mean absolute deviation
end
