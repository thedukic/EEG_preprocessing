function preproc_cleaning2(myPaths,id)
%
% Script for basic postprocessing of the preprocessed EEG data
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, January 2025
% =========================================================================

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject          = [];
subject.id       = id;
subject.task     = myPaths.task;
subject.group    = myPaths.group;
subject.visit    = myPaths.visit;
subject.preproc  = fullfile(myPaths.preproc, subject.id);
subject.clnfile1 = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'a.mat'];
subject.clnfile2 = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat'];

% Log time
t0 = datetime("now");
procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};

% Report
fprintf('\n');
disp('==================================================================');
fprintf('%s | %s | %s dataset | processing part 2\n',myPaths.group,subject.id,myPaths.task);
disp('==================================================================');
fprintf('\n');

% =========================================================================
% Basic post-preprocessing
% =========================================================================
% 1. Load cleaned data
fileName = fullfile(subject.preproc,subject.clnfile1);
if exist(fileName, 'file') == 2
    fprintf('%s: Loading the preprocessed data...\n',subject.id);
    load(fileName,'EEG');
else
    warning('%s: Cannot find the preprocessed data. Skipping...\n',subject.id);
    return;
end

% 2. Filter lowpass, must be done on continuous data
fprintf('\nLowpass filtering...\n');
EEG = do_filtering(EEG,'lowpass',cfg.flt);

% 3. IAF, should be done on continuous data?
fprintf('\nEstimating individual alpha frequency...\n');
EEG = check_restingIAF(EEG);

% 4. Epoch
fprintf('\nEpoching data...\n');
EEG = epoch_data(EEG,cfg.trg);

% 5. Common-average referening
EEG = do_reref(EEG,'aRegular');

% 6A.Traditional baseline correction
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    EEG = pop_rmbase(EEG,[],[]);
else
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0],[]);
end

% 6B. Regression based baseline correction method (recommended?)
% if strcmpi(myPaths.task,'MMN')
%     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0);
%
% elseif strcmpi(myPaths.task,'SART')
%     % SART wrt visual stimuli
%     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0);
%     % SART wrt response times
%     condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0);
%
% elseif strcmpi(myPaths.task,'MT')
%     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0);
%
% end
% if strcmpi(myPaths.task,'MMN') || strcmpi(myPaths.task,'SART')
%     EEG = correct_baseline(EEG,[(EEG.xmin)*1000 0],'Factor_1_Level_1',condLabel);
% else
%     error('Check this step for MT data.');
%     % How to consider MT2 and MT3 as one condition and MT5 as the other
%     % Make a copy and rename MT3 to MT2 in the .event struct?
%     % EEG = correct_baseline(EEG,[(EEG.xmin)*1000 0]); % if only 1 stimulus condition present
% end

% =========================================================================
% Data quaility checks
% =========================================================================
% Note the number of trials
NumberTrials = NaN(4,1);
NumberTrials(1) = size(EEG.data,3);

% =============================================
% A. EEGLAB-based rejection
% Any one of these functions can be commented out to ignore those artifacts when creating the mask
% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs

% ROIidx = 1:128; % Use only EEG electrodes!
% fprintf('\n================================\n');
% fprintf('Max. amplitude (>abs(%d uV))\n',cfg.epoch.rejectAmp);
% fprintf('================================\n');
% EEG = pop_eegthresh(EEG,1,ROIidx,-cfg.epoch.rejectAmp,cfg.epoch.rejectAmp,EEG.xmin,EEG.xmax,1,0);
%
% fprintf('\n================================\n');
% fprintf('Improbable data\n');
% fprintf('================================\n');
% EEG = pop_jointprob(EEG,1,ROIidx,cfg.epoch.singleChannelImprobableDataThreshold,cfg.epoch.allChannelImprobableDataThreshold,1,0);
%
% fprintf('\n================================\n');
% fprintf('Kurtosis\n');
% fprintf('================================\n');
% EEG = pop_rejkurt(EEG,1,ROIidx,cfg.epoch.singleChannelKurtosisThreshold,cfg.epoch.allChannelKurtosisThreshold,1,0);
%
% fprintf('\n================================\n');
% fprintf('Combining and rejecting\n');
% fprintf('================================\n');
% EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
% EEG = pop_rejepoch(EEG, EEG.reject.rejglobal, 0);

fprintf('\n================================\n');
fprintf('Max. amplitude (>abs(%d uV))\n',cfg.epoch.rejectAmp);
fprintf('================================\n');
% It could be that some people have very strong oscillations
% Example: Resting-state alpha oscillations
% We do not want to exclude these data
% -> bandstop filter 5-25 Hz
% -> EOG <5 Hz
% -> EMG >25 Hz
freqStop = [2 20]; % Hz

fprintf('Temporarily bandstop filtering [%d-%d Hz] the data.\n',freqStop);
fprintf('This prevents removal of trials with strong brain oscillations.\n');
eegchans = strcmp({EEG.chanlocs.type},'EEG');
dataTmp  = double(EEG.data(eegchans,:));
[b, a]   = butter(5, freqStop / (EEG.srate / 2), 'stop');
assert(isstable(b,a), 'Bandstop filter unstable.');
dataTmp  = filtfilt(b, a, dataTmp')';
dataTmp  = dataTmp - mean(dataTmp,1);
dataTmp  = dataTmp - mean(dataTmp,2);
dataTmp  = reshape(dataTmp, size(EEG.data(eegchans,:,:)));
badTrialTreshold = squeeze(any(abs(dataTmp) > cfg.epoch.rejectAmp, [1 2]));

% badTrialTmp = find(badTrialTreshold);
% figure; tiledlayout(1,2);
% mytopoplot(mean(dataTmp(:,:,badTrialTmp).^2,[2 3]),[],'Filtered',nexttile); colorbar;
% dataTmp = double(EEG.data(eegchans,:,:));
% mytopoplot(mean(dataTmp(:,:,badTrialTmp).^2,[2 3]),[],'Raw',nexttile); colorbar;
% disp(sum(badTrialTreshold));
%
% [NCHN,NPTS,NTRL]= size(dataTmp);
% for i = 1:NTRL
%     [psdspectra(:,:,i), freq] = pwelch(dataTmp(:,:,i)',NPTS,0,NPTS,EEG.srate);
% end
% figure; plot(freq, mean(psdspectra,3));

if any(badTrialTreshold)
    EEG = pop_rejepoch(EEG,badTrialTreshold,0);
else
    fprintf('No high voltage trials are found. \n');
end

% Note the number of trials
NumberTrials(2) = EEG.trials;

% =============================================
% B. EMG-slope-based rejection
fprintf('\n================================\n');
fprintf('EMG slopes\n');
fprintf('================================\n');

% Estimate slopes
slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
slopesChannelsxEpochs = slopesChannelsxEpochs > cfg.bch.muscleSlopeThreshold;

% tmp = any(slopesChannelsxEpochs,2);
% figure; mytopoplot(ones(1,128), tmp,'',nexttile);

% 1. Interpolate
% Interpolate trials that do not have a lot of electrodes contaminated by EMG
badElecsPerTrial = arrayfun(@(col) find(slopesChannelsxEpochs(:, col)), 1:size(slopesChannelsxEpochs, 2), 'UniformOutput', false);
eegchans = strcmp({EEG.chanlocs.type},'EEG');
chanLocs = EEG.chanlocs(eegchans);

% Max number of contaminated electrodes
fprintf('The maximum number of EMG-contaminated electrodes: %d\n', cfg.epoch.NinterpMax);
[data, report] = interpolate_epochs(EEG.data(eegchans,:,:),chanLocs,badElecsPerTrial,[],cfg.epoch.NinterpMax);

% Return the data to the struct
EEG.data(eegchans,:,:) = data;

% These are too noisy to be saved
badTrialMuscle = false(size(badElecsPerTrial));
badTrialMuscle(report.listNotFixed) = true;

% 2. Remove all
% BadEpochs = sum(slopesChannelsxEpochs, 1);
% badTrialMuscleTreshhold = 0;
% badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;
%
% while sum(badTrialMuscle) > 0.5*EEG.trials
%     warning('Increasing the minimum number of allowed EMG-contaminated channels: %d -> %d!',badTrialMuscleTreshhold,badTrialMuscleTreshhold+1);
%     badTrialMuscleTreshhold = badTrialMuscleTreshhold+1;
%     badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;
% end

if any(badTrialMuscle)
    EEG = pop_rejepoch(EEG,badTrialMuscle,0);
else
    fprintf('No EMG-contaminated trials are found.\n');
end

% Log
EEG.ALSUTRECHT.epochRejections.InterpTrialInfo = badElecsPerTrial;
EEG.ALSUTRECHT.epochRejections.InterpReport = report;

% Note the number of trials
NumberTrials(3) = EEG.trials;

% =============================================
% % C. Detection using variance and the G-ESD method
% fprintf('\n================================\n');
% fprintf('Variance and the G-ESD method\n');
% fprintf('================================\n');
% EEG = detect_badepochs(EEG);

% Note the number of trials
NumberTrials(4) = EEG.trials;

% =========================================================================
fprintf('\n================================\n');
fprintf('Final checks\n');
fprintf('================================\n');

% % IAF
% fprintf('Estimating individual alpha frequency...\n');
% EEG = check_restingIAF(EEG);

% Median voltage shift
fprintf('Estimating voltage range...\n');
voltageShiftWithinEpoch = median(range(EEG.data,2),3);
EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal = voltageShiftWithinEpoch;

% Channel correlation matrix
fprintf('Estimating channel correlation matrix...\n');
EEG = check_channelcov(EEG);

% % Empirical frequency boundaries
% fprintf('Estimating empirical frequency boundaries...\n');
% EEG = estimate_gedBounds(EEG);

% EMG leftovers
fprintf('Checking EMG leftovers...\n');
slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
slopesChannelsxEpochs(slopesChannelsxEpochs < cfg.bch.muscleSlopeThreshold) = NaN;
slopesChannelsxEpochs = slopesChannelsxEpochs - cfg.bch.muscleSlopeThreshold;
BadEpochs = sum(slopesChannelsxEpochs,1,'omitnan');
muscleLeftover = mean(BadEpochs > 0);
fprintf('EMG leftovers: %1.2f\n',muscleLeftover);

% Log
EEG.ALSUTRECHT.epochRejections.initialEpochs       = NumberTrials(1);
EEG.ALSUTRECHT.epochRejections.afterEEGLABEpochs   = NumberTrials(2);
EEG.ALSUTRECHT.epochRejections.afterEMGSlopeEpochs = NumberTrials(3);
EEG.ALSUTRECHT.epochRejections.remainingEpochs     = NumberTrials(4);
EEG.ALSUTRECHT.epochRejections.interpEpochs        = length(report.listFixed);
EEG.ALSUTRECHT.epochRejections.proportionOfEpochsRejected = (EEG.ALSUTRECHT.epochRejections.initialEpochs-EEG.ALSUTRECHT.epochRejections.remainingEpochs)/EEG.ALSUTRECHT.epochRejections.initialEpochs;
EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal = voltageShiftWithinEpoch;

% EEG.ALSUTRECHT.epochRejections.badTrialMuscleTreshhold2 = badTrialMuscleTreshhold;
EEG.ALSUTRECHT.epochRejections.muscle2 = muscleLeftover;
EEG.ALSUTRECHT.leftovers.muscle2 = EEG.ALSUTRECHT.epochRejections.muscle2;

% =========================================================================
% Report
fprintf('\n================================\n');
fprintf('Final overview\n');
fprintf('================================\n');
% Nremoved = [NumberTrials(1)-NumberTrials(2), NumberTrials(2)-NumberTrials(3)];
Nremoved = abs(diff(NumberTrials));
fprintf('Inital trials:    %d\n',NumberTrials(1));
fprintf('Removed trials:   %d + %d + %d = %d\n',Nremoved,sum(Nremoved));
fprintf('Remaining trials: %d\n',NumberTrials(end));

% =========================================================================
% Estimate the spectra
[psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG,'preproc2');

% Plots differe per task
if strcmpi(myPaths.task,'MMN') || strcmpi(myPaths.task,'SART')
    maskCond = [EEG.event.edftype];
    if strcmpi(myPaths.task,'SART')
        condTrig = [3 6];
        maskTrig = maskCond==condTrig(1) | maskCond==condTrig(2);
        minClim  = [-8 8];
        NTLS     = 3;
    elseif strcmpi(myPaths.task,'MMN')
        condTrig = [17 12];
        maskTrig = maskCond==condTrig(1) | maskCond==condTrig(2);
        minClim  = [-3 3];
        NTLS     = 4;
    end

    NTRG = length(condTrig);
    maskCond = maskCond(maskTrig);
    assert(EEG.trials==sum(maskTrig));
    % condTrig = unique(maskCond);

    % Select only EEG
    dataCmap = brewermap(128,'PRGn');
    % dataCmap = brewermap(128,'BrBG');

    fh = figure;
    th = tiledlayout(1,NTLS);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    dataTmp1 = mean(EEG.data(chaneeg,:,maskCond==condTrig(1)),3);
    dataTmp2 = mean(EEG.data(chaneeg,:,maskCond==condTrig(2)),3);

    % ERP 1&2
    for i = 1:NTRG
        if i == 1
            dataTmp = dataTmp1;
        else
            dataTmp = dataTmp2;
        end

        % dataTmp = mean(EEG.data(chaneeg,:,maskCond==condTrig(i)),3);
        NumberTrials = sum(maskCond==condTrig(i));

        dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
        dataClim(1) = floor(dataClim(1));
        dataClim(2) = ceil(dataClim(2));

        dataClim(1) = min(dataClim(1),minClim(1));
        dataClim(2) = max(dataClim(2),minClim(1));

        th = nexttile;
        hold on; box off;
        plot([0 0],dataClim,'Color',0.7*ones(1,3),'LineWidth',1,'HandleVisibility','off');
        plot(EEG.times([1 end]),[0 0],'Color',0.7*ones(1,3),'LineWidth',1,'HandleVisibility','off');
        plot(EEG.times,dataTmp,'LineWidth',1.1);

        colororder(th,dataCmap); xlim(EEG.times([1 end])); ylim(dataClim);
        title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, trig ' num2str(condTrig(i)) ', N = ' num2str(NumberTrials) ' trials']);
        pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');
    end

    % ERP difference
    dataTmp = dataTmp1 - dataTmp2;

    dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));

    dataClim(1) = min(dataClim(1),minClim(1));
    dataClim(2) = max(dataClim(2),minClim(1));

    th = nexttile;
    hold on; box off;
    plot([0 0],dataClim,'Color',0.7*ones(1,3),'LineWidth',1,'HandleVisibility','off');
    plot(EEG.times([1 end]),[0 0],'Color',0.7*ones(1,3),'LineWidth',1,'HandleVisibility','off');
    plot(EEG.times,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap); xlim(EEG.times([1 end])); ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, ERP difference (' num2str(condTrig(1)) '-' num2str(condTrig(2)) ')']);
    pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');

    % Add MMN topolot
    if strcmpi(myPaths.task,'MMN')
        X  = mean(dataTmp(chaneeg,EEG.times>150 & EEG.times<300),2);
        X0 = mean(dataTmp(chaneeg,EEG.times<0),2);
        [bl, al] = butter(2,20/(EEG.srate/2),'low'); assert(isstable(bl,al));
        XX = filtfilt(bl,al, X-X0);
        mytopoplot(XX,[],'Lowpass filtered MMN (<20 Hz, 150-300 ms)',nexttile,0.5*[-1 1]);
    end

    plotX=30; plotY=8;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[subject.id '_erp_final']),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra,2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

elseif strcmpi(myPaths.task,'RS')
    % =============================
    % Plot 1
    % =============================
    dataCmap = brewermap(sum(chaneeg),'BrBG');

    fh = figure;
    th = tiledlayout(1,2);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    th = nexttile;
    hold on; box off;
    plot(freq,psdspectra,'LineWidth',1.1);

    colororder(th,dataCmap); xlim([freq(1), 60]); % ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('Power (a.u.)');

    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra);
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freq,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap); xlim(freq([1 end])); ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('log_{10}(Power) (a.u.)');

    plotX=20; plotY=8;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[subject.id '_pspectra_final']),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra,2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

elseif strcmpi(myPaths.task,'MT')
    dataCmap1 = brewermap(sum(chaneeg),'BrBG');
    dataCmap2 = brewermap(sum(chanemg),'PRGn');

    fh = figure;
    th = tiledlayout(2,2);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    % EEG
    th = nexttile;
    hold on; box off;
    plot(freq,psdspectra(:,chaneeg),'LineWidth',1.1);

    colororder(th,dataCmap1); xlim([freq(1), 60]); % ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('Power (a.u.)');

    % EMG
    th = nexttile;
    hold on; box off;
    plot(freq,psdspectra(:,chanemg),'LineWidth',1.1);

    colororder(th,dataCmap2); xlim(freq([1 end])); % ylim(dataClim);
    title([subject.id ', ' num2str(sum(chanemg)) ' EMG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('Power (a.u.)');

    % log(EEG)
    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra(:,chaneeg));
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freq,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap1); xlim(freq([1 end])); ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('log_{10}(Power) (a.u.)');

    % log(EMG)
    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra(:,chanemg));
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freq,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap2); xlim(freq([1 end])); ylim(dataClim);
    title([subject.id ', ' num2str(sum(chanemg)) ' EMG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('log_{10}(Power) (a.u.)');

    plotX=20; plotY=14;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[subject.id '_pspectra_final']),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra(:,chaneeg),2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

end

% =========================================================================
% Remove IC signals
EEG.icaact = [];

% Save
fprintf('\n%s: Saving the preprocessed data (part 2)...\n',subject.id);
save(fullfile(subject.preproc,subject.clnfile2),'EEG','cfg','procTimeTags');

% Report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf('Finished: %s\n',t1);
fprintf('Running time: %d min.\n\n',dd);

end