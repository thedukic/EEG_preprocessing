function checkReport = preproc_cleaning2(myPaths,id)
%
% Script for basic postprocessing of the preprocessed EEG data
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, August 2024
% =========================================================================
% TODO
% 1. Lowpass filtering for ERPs (<30-40 Hz) and epoching to be done here
% 2.
%

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject          = [];
subject.id       = id;
subject.task     = myPaths.task;
subject.group    = myPaths.group;
subject.visit    = myPaths.visit;
subject.preproc  = fullfile(myPaths.preproc, subject.id);
subject.clnfile  = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' cfg.rnum '.mat'];
subject.clnfile1 = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' cfg.rnum 'b.mat'];

% Log time
t0 = datetime("now");
procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};

% Report
fprintf('\n');
disp('==================================================================');
fprintf('%s | %s | %s dataset | processing part 2\n',myPaths.group,subject.id,myPaths.task);
disp('==================================================================');
fprintf('\n');

% 1. Load cleaned data
fileName = fullfile(subject.preproc,subject.clnfile);
if exist(fileName, 'file') == 2
    fprintf('%s: Loading the preprocessed data...\n',subject.id);
    load(fileName,'EEG');
else
    warning('%s: Cannot find the preprocessed data. Skipping...\n',subject.id);
    checkReport = [];
    return;
end


% 2. Common-average referening and baseline correction
EEG = do_reref(EEG,'aRegular');

% 2A.Traditional baseline correction
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    EEG = pop_rmbase(EEG,[],[]);
else
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0],[]);
end

% 2B. Regression based baseline correction method (recommended)
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

% 3. Reject epochs that are likely still noisy
NumberTrials = NaN(3,1);
NumberTrials(1) = size(EEG.data,3);

% Any one of these functions can be commented out to ignore those artifacts when creating the mask
% This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs
% A. EEGLAB-based rejection
ROIidx = 1:128;
fprintf('\n================================\n');
fprintf('Max. amplitude\n');
fprintf('================================\n');
EEG = pop_eegthresh(EEG,1,ROIidx,-cfg.epoch.rejectAmp,cfg.epoch.rejectAmp,EEG.xmin,EEG.xmax,1,0);
fprintf('\n================================\n');
fprintf('Improbable data\n');
fprintf('================================\n');
EEG = pop_jointprob(EEG,1,ROIidx,cfg.epoch.singleChannelImprobableDataThreshold,cfg.epoch.allChannelImprobableDataThreshold,1,0);
fprintf('\n================================\n');
fprintf('Kurtosis\n');
fprintf('================================\n');
EEG = pop_rejkurt(EEG,1,ROIidx,cfg.epoch.singleChannelKurtosisThreshold,cfg.epoch.allChannelKurtosisThreshold,1,0);

fprintf('\n================================\n');
fprintf('Combining and rejecting\n');
fprintf('================================\n');
EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
EEG = pop_rejepoch(EEG, EEG.reject.rejglobal, 0);

% B. EMG-slope-based rejection
fprintf('\n================================\n');
fprintf('EMG slopes\n');
fprintf('================================\n');
slopesChannelsxEpochs = dected_emg(EEG);
slopesChannelsxEpochs = slopesChannelsxEpochs > cfg.bch.muscleSlopeThreshold;

% How many channels are affected in each epoch
BadEpochs = sum(slopesChannelsxEpochs, 1);
badTrialMuscleTreshhold = 0;
badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;

NumberTrials(2) = EEG.trials;
while sum(badTrialMuscle) > 0.5*EEG.trials
    warning('Increasing the minimum number of allowed EMG-contaminated channels: %d -> %d!',badTrialMuscleTreshhold,badTrialMuscleTreshhold+1);
    badTrialMuscleTreshhold = badTrialMuscleTreshhold+1;
    badTrialMuscle = BadEpochs>badTrialMuscleTreshhold;
end

EEG = pop_rejepoch(EEG,badTrialMuscle,0);
NumberTrials(3) = EEG.trials;

% 4. Check the median voltahe shift
% The following measures the median size of the voltage shift in each epoch for each participant.
% This is used in 'report_final' to provide advice on which participants to manually check
voltageShiftWithinEpoch = range(EEG.data(:,:,:),2);

% 5. Log
EEG.ALSUTRECHT.epochRejections.initialEpochs   = NumberTrials(1);
EEG.ALSUTRECHT.epochRejections.round1Epochs    = NumberTrials(2);
EEG.ALSUTRECHT.epochRejections.remainingEpochs = NumberTrials(3);
EEG.ALSUTRECHT.epochRejections.proportionOfEpochsRejected = (EEG.ALSUTRECHT.epochRejections.initialEpochs-EEG.ALSUTRECHT.epochRejections.remainingEpochs)/EEG.ALSUTRECHT.epochRejections.initialEpochs;

EEG.ALSUTRECHT.issues_to_check.MedianvoltageshiftwithinepochFinal = median(voltageShiftWithinEpoch,3);
EEG.ALSUTRECHT.issues_to_check.NumberTrials3 = EEG.ALSUTRECHT.epochRejections.round1Epochs;
EEG.ALSUTRECHT.issues_to_check.NumberTrials4 = EEG.ALSUTRECHT.epochRejections.remainingEpochs;

% Report
Nremoved = [NumberTrials(1)-NumberTrials(2), NumberTrials(2)-NumberTrials(3)];
fprintf('\n%s: Removed trials %d + %d = %d\n',subject.id,Nremoved,sum(Nremoved));
fprintf('%s: Remaining trials %d\n',subject.id,NumberTrials(3));

% 6. Plot
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
    chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
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
    dataTmp = dataTmp1-dataTmp2;

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

elseif strcmpi(myPaths.task,'RS')
    % Select only EEG
    chaneeg = strcmp({EEG.chanlocs.type},'EEG');
    dataeeg = EEG.data(chaneeg,:,:);

    % Compute power spectra
    [NCHN,NPTS,NTRL]= size(dataeeg);
    psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

    for i = 1:NTRL
        [psdspectra(:,:,i), freq] = pwelch(dataeeg(:,:,i)',NPTS,0,NPTS,EEG.srate);
    end

    % Average the spectra
    psdspectra = mean(psdspectra,3);

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

    dataTmp = log10(psdspectra);
    dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
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

elseif strcmpi(myPaths.task,'MT')
    % Select only EEG/EMG
    chaneeg = strcmp({EEG.chanlocs.type},'EEG');
    chanemg = strcmp({EEG.chanlocs.type},'EMG');
    dataeeg = EEG.data;

    % Compute power spectra
    [NCHN,NPTS,NTRL]= size(dataeeg);
    psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

    for i = 1:NTRL
        [psdspectra(:,:,i), freq] = pwelch(dataeeg(:,:,i)',NPTS,0,NPTS,EEG.srate);
    end

    % Average the spectra
    psdspectra = mean(psdspectra,3);

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

    dataTmp = log10(psdspectra(:,chaneeg));
    dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freq,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap1); xlim(freq([1 end])); ylim(dataClim);
    title([subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequnecy (Hz)'); ylabel('log_{10}(Power) (a.u.)');

    % log(EMG)
    th = nexttile;
    hold on; box off;

    dataTmp = log10(psdspectra(:,chanemg));
    dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
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
end

% 7. Save and return
fprintf('%s: Saving the preprocessed data (part 2)...\n',subject.id);
checkReport = EEG.ALSUTRECHT.issues_to_check;
checkReport = rmfield(checkReport,{'Medianvoltageshiftwithinepoch','MedianvoltageshiftwithinepochFinal'});

EEG.icaact = [];
preprocReport = EEG.ALSUTRECHT;
save(fullfile(subject.preproc,subject.clnfile1),'EEG','preprocReport','cfg','procTimeTags');

end