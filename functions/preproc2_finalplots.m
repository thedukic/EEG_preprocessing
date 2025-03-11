function preproc2_finalplots(EEG,subject,NumberTrials,thisTask,thisTag)

fprintf('\n================================\n');
fprintf('Final reports\n');
fprintf('================================\n');

% Report
% Nremoved = [NumberTrials(1)-NumberTrials(2), NumberTrials(2)-NumberTrials(3)];
Nremoved = abs(diff(NumberTrials));
fprintf('Inital trials:    %d\n',NumberTrials(1));
fprintf('Removed trials:   %d + %d + %d = %d\n',Nremoved,sum(Nremoved));
fprintf('Remaining trials: %d\n',NumberTrials(end));

% Estimate power sepctra
[psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG,'preproc2');
freqlog = log10(freq);

% Freq labels
freqTicks     = [0:6 8 10 15 20 30 50 100];
freqTicksLog  = log10(freqTicks);
freqTicksCell = arrayfun(@num2str, freqTicks, 'UniformOutput', false);

% Plots differe per task
if strcmpi(thisTask,'MMN') || strcmpi(thisTask,'SART')
    maskCond = [EEG.event.edftype];
    if strcmpi(thisTask,'SART')
        if strcmpi(EEG.ALSUTRECHT.SART.type,'StimulusLocked')
            condTrig = [3 6];
            maskTrig = maskCond==condTrig(1) | maskCond==condTrig(2);
            NTLS     = 3;
        elseif strcmpi(EEG.ALSUTRECHT.SART.type,'ResponseLocked')
            condTrig = 1;
            % maskTrig = maskCond==condTrig;
            NTLS     = 1;
        end
        minClim = [-8 8];

    elseif strcmpi(thisTask,'MMN')
        condTrig = [17 12];
        maskTrig = maskCond==condTrig(1) | maskCond==condTrig(2);
        minClim  = [-3 3];
        NTLS     = 4;
    end

    % Number of conditions
    NTRG = length(condTrig);

    if NTRG == 2
        maskCond = maskCond(maskTrig);
        assert(EEG.trials == sum(maskTrig));
        % condTrig = unique(maskCond);

        % ERP 1/2
        dataTmp1 = mean(EEG.data(chaneeg,:,maskCond==condTrig(1)),3);
        dataTmp2 = mean(EEG.data(chaneeg,:,maskCond==condTrig(2)),3);

        NumberTrials = NaN(NTRG,1);
        NumberTrials(1) = sum(maskCond == condTrig(1));
        NumberTrials(2) = sum(maskCond == condTrig(2));

    elseif NTRG == 1
        % ERP 1
        dataTmp1 = mean(EEG.data(chaneeg,:,:),3);
        NumberTrials = EEG.trials;

    end

    % Select only EEG
    dataCmap = brewermap(128,'PRGn');
    % dataCmap = brewermap(128,'BrBG');

    fh = figure;
    th = tiledlayout(1,NTLS);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    % ERP 1&2
    for i = 1:NTRG
        if i == 1
            dataTmp = dataTmp1;
        else
            dataTmp = dataTmp2;
        end

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
        title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, trig ' num2str(condTrig(i)) ', N = ' num2str(NumberTrials(i)) ' trials']);
        pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');
    end

    % ERP difference
    if NTRG == 2
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
        title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, ERP difference (' num2str(condTrig(1)) '-' num2str(condTrig(2)) ')']);
        pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');
    end

    % Add MMN topolot
    if strcmpi(thisTask,'MMN')
        X  = mean(dataTmp(chaneeg,EEG.times>150 & EEG.times<300),2);
        X0 = mean(dataTmp(chaneeg,EEG.times<0),2);
        [bl, al] = butter(2,20/(EEG.srate/2),'low'); assert(isstable(bl,al));
        XX = filtfilt(bl,al, X-X0);
        mytopoplot(XX,[],'Lowpass filtered MMN (<20 Hz, 150-300 ms)',nexttile,0.5*[-1 1]);
    end

    plotX=30; plotY=8;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[EEG.ALSUTRECHT.subject.id '_erp_final_' num2str(thisTag)]),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra,2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

elseif strcmpi(thisTask,'RS')
    % =============================
    % Plot 1
    % =============================
    dataCmap = brewermap(sum(chaneeg),'BrBG');

    fh = figure;
    th = tiledlayout(1,2);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    % Regular plot
    th = nexttile;
    hold on; box off;
    plot(freq,psdspectra,'LineWidth',1.1);

    colororder(th,dataCmap); xlim([freq(1), 60]); % ylim(dataClim);
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');

    % Log-log plot
    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra);
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freqlog,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap); xlim(freqlog([1 end])); ylim(dataClim);
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
    xticks(freqTicksLog); xticklabels(freqTicksCell);

    plotX=20; plotY=8;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[EEG.ALSUTRECHT.subject.id '_pspectra_final_' num2str(thisTag)]),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra,2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

elseif strcmpi(thisTask,'MT')
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
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');

    % EMG
    th = nexttile;
    hold on; box off;
    plot(freq,psdspectra(:,chanemg),'LineWidth',1.1);

    colororder(th,dataCmap2); xlim(freq([1 end])); % ylim(dataClim);
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chanemg)) ' EMG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');

    % log(EEG)
    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra(:,chaneeg));
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freqlog,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap1); xlim(freqlog([1 end])); ylim(dataClim);
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chaneeg)) ' EEG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
    xticks(freqTicksLog); xticklabels(freqTicksCell);

    % log(EMG)
    th = nexttile;
    hold on; box off;

    dataTmp     = log10(psdspectra(:,chanemg));
    dataClim    = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
    dataClim(1) = floor(dataClim(1));
    dataClim(2) = ceil(dataClim(2));
    plot(freqlog,dataTmp,'LineWidth',1.1);

    colororder(th,dataCmap2); xlim(freqlog([1 end])); ylim(dataClim);
    title([EEG.ALSUTRECHT.subject.id ', ' num2str(sum(chanemg)) ' EMG, N = ' num2str(NumberTrials(3)) ' trials']);
    pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
    xticks(freqTicksLog); xticklabels(freqTicksCell);

    plotX=20; plotY=14;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(subject.preproc,[EEG.ALSUTRECHT.subject.id '_pspectra_final_' num2str(thisTag)]),'-dtiff','-r300');
    close(fh);

    % =============================
    % Plot 2 (Freq bands and IAF)
    % =============================
    plot_iaf(EEG.ALSUTRECHT.pSpec.sums.paf,mean(psdspectra(:,chaneeg),2),freq,EEG.ALSUTRECHT.pSpec.sums.muSpec,EEG.ALSUTRECHT.pSpec.sums.freq,subject);

end

% fprintf('Done!\n');

end