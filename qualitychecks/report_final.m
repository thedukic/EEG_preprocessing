function report_final(myPaths,subjects)
%
% Script for reporting on EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, January 2025
% =========================================================================
% TODO
% 1.
% 2.
%

fprintf('\n');
disp('==================================================================');
fprintf('Generating the final report for %s. This may take a bit...\n',myPaths.group);
disp('==================================================================');
fprintf('\n');

% =========================================================================
NSUB = length(subjects);
chanlocs = readlocs('biosemi128_eeglab.ced');
chanlbls = {chanlocs.labels};
maskelec = zeros(length(chanlocs),NSUB);

% Report folder
myPaths.reports = fullfile(myPaths.preproc,'reports');
if exist(myPaths.reports,'dir')~=7, mkdir(myPaths.reports); end

% =========================================================================
N = NaN(NSUB,5);
NCHN = length(chanlbls);
Medianvoltageshift = NaN(NCHN,NSUB);
CorrelationMatrices = NaN(NCHN+2,NCHN+2,NSUB);

fprintf('Loading %d datasets. Wait...\n',NSUB);
for i = 1:NSUB
    load(fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat']),'EEG');

    % Record warnings about potential issues
    EEG = report_issues(EEG);

    % Record warnings for all participants in a single table
    if i == 1
        tableIssues = struct2table(EEG.ALSUTRECHT.issues_to_check,'AsArray',true);
        tableIssues(2:NSUB,:) = cell2table([repmat({''},NSUB-1,1), repmat({NaN},NSUB-1,length(tableIssues.Properties.VariableNames)-1)], 'VariableNames', tableIssues.Properties.VariableNames);
    else
        tableIssues(i,:) = struct2table(EEG.ALSUTRECHT.issues_to_check,'AsArray',true);
    end

    % Bad electrodes
    maskelec(:,i) = double(ismember(chanlbls,EEG.ALSUTRECHT.badchaninfo.badElectrodes));
    N(i,1) = sum(maskelec(:,i));

    % Epochs: Total possible
    % N(i,2) = sum([EEG.ALSUTRECHT.eventinfo{:,3}])*4/2; % RS
    N(i,2) = EEG.ALSUTRECHT.issues_to_check.NumberTrials1;
    % Epochs: Left after preproc1
    N(i,3) = EEG.ALSUTRECHT.issues_to_check.NumberTrials2;
    % Epochs: Left after preproc2
    N(i,4) = EEG.ALSUTRECHT.issues_to_check.NumberTrials3;

    % Leftover EMG
    % N(i,5) = EEG.ALSUTRECHT.leftovers.muscle1;  % after proc1
    N(i,5) = EEG.ALSUTRECHT.leftovers.muscle2;    % after proc2

    % Median voltage range
    Medianvoltageshift(:,i) = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

    % Correlation matrix
    CorrelationMatrices(:,:,i) = EEG.ALSUTRECHT.chanCorr;
end

% =========================================================================
% Plot
fh = figure;
th = tiledlayout(4,1);
th.TileSpacing = 'compact'; th.Padding = 'compact';

% Plot 1: Which channels are usually interpoalted
myCmap1 = brewermap(128,'RdPu');
myCmap2 = brewermap(3,'YlOrRd');
myCmap3 = brewermap(12,'Paired');

nexttile(1);
topoplot(mean(maskelec,2),chanlocs,'maplimits',[0 0.25],'headrad',0.5,'colormap',myCmap1,'whitebk','on','electrodes','on','style','map','shading','interp');
axis tight; title([myPaths.group ', N = ' num2str(NSUB)]);
hcb = colorbar;
hcb.Title.String = "%";

% Plot 2: Interpolated channels per person
nexttile; hold on;
plot([0 NSUB+1],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,1)/128,'FaceColor',myCmap1(end/2,:)); ylim([0 1]); box on;

ylabel('Interp. chans (%)');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);

% Plot 3: EMG lefovers from preprocessing
nexttile; hold on;
plot([0 NSUB+1],[0.25 0.25],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,5),'FaceColor',myCmap3(12,:)); ylim([0 1]); box on;

xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);
ylabel('EMG leftover (%)');

% Plot 4: Number of trials overview
nexttile;
Ntmp = [N(:,4), N(:,3)-N(:,4), N(:,2)-N(:,3)];
bh = bar(Ntmp,'stacked');
bh(1).FaceColor = myCmap2(3,:);
bh(2).FaceColor = myCmap2(2,:);
bh(3).FaceColor = myCmap2(1,:);

ylim([0 round(max(sum(Ntmp,2))*1.1)]);
ylabel('# trials');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);
% xlabel('Participant');

% Save
plotX=25; plotY=22;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myPaths.reports,['Summary1_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r400');

% =========================================================================
% The following checks for participants who show outlying values for the median voltage shift within each epoch:
% The following detects outlier files in the median amount of their max-min
% voltage shift within an epoch, after adjusting for the fact that the data
% across all participants is likely to be positively skewed with a log transform.
MedianvoltageshiftwithinepochLogged = log10(Medianvoltageshift);
InterQuartileRange = iqr(MedianvoltageshiftwithinepochLogged,2);
Upper25 = prctile(MedianvoltageshiftwithinepochLogged,75,2);
Lower25 = prctile(MedianvoltageshiftwithinepochLogged,25,2);

% 75th% and 25th% +/- (1.5 x IQR) is the recommended outlier detection method,
% so this is used to recommend which participants to manually check
% However, I find this can be a bit too sensitive upon manual inspection,
% and that 1.75, 2, or even 2.5 can be a better threshold.
LowerBound = NaN(128,1);
UpperBound = NaN(128,1);
for i = 1:size(MedianvoltageshiftwithinepochLogged,1)
    LowerBound(i) = Lower25(i) - (2 * InterQuartileRange(i));
    UpperBound(i) = Upper25(i) + (2 * InterQuartileRange(i));
end

VoltageShiftsTooLow = MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooLow = VoltageShiftsTooLow - LowerBound;
VoltageShiftsTooLow(VoltageShiftsTooLow > 0) = 0;
CumulativeSeverityOfAmplitudesBelowThreshold = sum(VoltageShiftsTooLow,1)';

VoltageShiftsTooHigh = MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooHigh = VoltageShiftsTooHigh - UpperBound;
VoltageShiftsTooHigh(VoltageShiftsTooHigh < 0) = 0;
CumulativeSeverityOfAmplitudesAboveThreshold = sum(VoltageShiftsTooHigh,1)';

% 2. Plot:
fh = figure; hold on;
plot(LowerBound,'LineWidth',1.5,'Color','k');
plot(UpperBound,'LineWidth',1.5,'Color','k');
plot(MedianvoltageshiftwithinepochLogged,'LineWidth',1.2);

ylabel('Median voltage shift (uV)');
xticks(1:128); xticklabels({chanlocs.labels}); xlim([1 128]);
legend('LowerBound', 'UpperBound');
set(gca,'ColorOrder',[0 0 0; 0 0 0; brewermap(NSUB,'BuGn')]);

plotX=35; plotY=15;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myPaths.reports,['Summary2_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r400');

% 3. Visualise
maskSubj = find(CumulativeSeverityOfAmplitudesAboveThreshold > 0);

if ~isempty(maskSubj)
    NSUBhv = length(maskSubj);

    fh = figure;
    th = tiledlayout(1,NSUBhv);
    th.TileSpacing = 'compact'; th.Padding = 'compact';
    title(th,{'Participants with high median voltage','Sometimes may be just very strong alpha waves'});

    MedianValue = prctile(MedianvoltageshiftwithinepochLogged,50,"all");
    myClim = [MedianValue, mean(UpperBound)];
    % myClim = [mean(LowerBound), mean(UpperBound)];

    for i = 1:NSUBhv
        maskChan = VoltageShiftsTooHigh(:,maskSubj(i))>0;
        str = strjoin(chanlbls(maskChan),', ');
        fprintf('%s: %s\n',subjects{maskSubj(i)},str);

        % [min(MedianvoltageshiftwithinepochLogged(:,maskSubj(i))), max(MedianvoltageshiftwithinepochLogged(:,maskSubj(i)))]
        nexttile;
        topoplot(MedianvoltageshiftwithinepochLogged(:,maskSubj(i)),chanlocs,'maplimits',myClim,'headrad', 0.5,'colormap',myCmap1,'whitebk','on','electrodes','on','style','map','emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(maskChan),'o','k',4,1},'shading','interp');
        title(subjects{maskSubj(i)});
    end

    plotX=25; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(myPaths.reports,['Summary3_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r400');

else
    fprintf('None of the participants in this group have an unusually high voltage.\n');
end

% OutlierParticipantsToManuallyCheck   = table(Participant_IDs', CumulativeSeverityOfAmplitudesBelowThreshold,CumulativeSeverityOfAmplitudesAboveThreshold);
% LoggedMedianVoltageShiftAcrossEpochs = array2table(MedianvoltageshiftwithinepochLogged);

% =========================================================================
% CorrelationMatricesMean = mean(CorrelationMatrices,3);
%
% NCHN = size(CorrelationMatricesMean,1);
% CorrelationMatricesZ = 0*CorrelationMatrices;
%
% for i = 1:NCHN
%     for j = 1:NCHN
%         CorrelationMatricesZ(i,j,:) = zscore(squeeze(CorrelationMatrices(i,j,:)));
%     end
% end
%
% maskDiff = squeeze(any(abs(CorrelationMatricesZ)>6,[1 2]));
% maskDiff = find(maskDiff);
%
% maskDiff = [1; 2; 3; maskDiff]
%
% NSUBdv = length(maskDiff);
%
% fh = figure;
% th = tiledlayout(1,NSUBdv+1);
% th.TileSpacing = 'compact'; th.Padding = 'compact';
%
% for i = 1:NSUBdv+1
%     nexttile;
%     if i == 1
%         imagesc(CorrelationMatricesMean);
%         title('Group average');
%     else
%         imagesc(CorrelationMatricesZ(:,:,maskDiff(i-1)));
%         title(subjects{maskDiff(i-1)});
%     end
%     axis square;
% end

% =========================================================================
% Empirical estimates
cutoffSpread = 0.04;
cutoffPeak   = 5;

cnt = 0;
psdSpreadsAll = NaN(NSUB,1);
psdPeaksAll   = NaN(NSUB,1);

fprintf('Loading %d datasets... Wait...\n',NSUB);
for i = 1:NSUB
    load(fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat']),'EEG');

    % Estimate the spectra
    [psdspectra, freq, chaneeg] = estimate_power(EEG,'freport');

    % =================================================================
    % 1. Calculate measure of power spread across EEG electrodes
    maskFreq  = freq>40 & freq<70;
    % psdSpread = psdspectra ./ sum(psdspectra,1);
    psdSpread = psdspectra;
    psdSpreadsAll(i) = mean(std(psdSpread(maskFreq,:),0,2));
    % psdSpreadsAll(i) = std(mean(psdspectra(mask,:),1));

    % =================================================================
    % 2. Find the number of peaks in the average EEG spectrum
    psdPeak  = mean(psdspectra,2);
    psdPeak  = psdPeak ./ sum(psdPeak,1);
    maskFreq = freq < 40;
    [qrspeaks, locs] = findpeaks(psdPeak(maskFreq),freq(maskFreq),'MinPeakProminence',0.0005);
    psdPeaksAll(i) = length(locs);

    % =================================================================
    % figure; plot(std(psdspectra(mask,:),0,2))
    % figure; histogram(std(psdspectra(mask,:),0,2))

    % Plot and output if it exceeds the cut-off
    if psdSpreadsAll(i) > cutoffSpread || psdPeaksAll(i) > cutoffPeak
        fprintf('%s : Spread = (%.2f), Peaks = %d.\n', subjects{i}, psdSpreadsAll(i), psdPeaksAll(i));

        cnt = cnt+1;
        if cnt == 1
            dataCmap = brewermap(sum(chaneeg),'BrBG');

            fh = figure;
            th = tiledlayout("flow");
            th.TileSpacing = 'compact'; th.Padding = 'compact';
        end

        th = nexttile;
        hold on; box off;

        % Plot log power
        psdspectra  = log10(psdspectra);
        % dataClim    = 1.1*[min(psdspectra(:)), max(psdspectra(:))];
        % dataClim(1) = floor(dataClim(1));
        % dataClim(2) = ceil(dataClim(2));
        plot(freq,psdspectra,'LineWidth',1.1);
        colororder(th,dataCmap);

        for j = 1:psdPeaksAll(i)
            plot(locs(j) * ones(1,2),[-4 3],'LineWidth',1,'Color',0.2*ones(1,3));
            % scatter(locs,psdspectra(any(freq==locs',2)),'filled','*','MarkerFaceColor',0.2*ones(1,3));
        end

        xlim([freq(1), 70]); ylim([-4 3]); % ylim(dataClim);
        title({[subjects{i},', group: ',myPaths.group, ', visit: ', myPaths.visit, ', task: ', myPaths.task], ['Spread: ', num2str(psdSpreadsAll(i), '%.2f') ', Peaks: ', num2str(psdPeaksAll(i), '%d') ]});
        pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('log_{10}(Power) (a.u.)');
    end
end

if cnt > 0
    plotX=25; plotY=25;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(myPaths.reports,['Summary4_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r300');
    % close(fh);
end

% Plot the distributions of the power spectra characteristics
% -> Power spectra spread, power spectra peaks
fh = figure;
th = tiledlayout(1,2);
th.TileSpacing = 'compact'; th.Padding = 'compact';
nexttile;
hb = histogram(psdSpreadsAll);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power spectra spread');
nexttile;
hb = histogram(psdPeaksAll);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power spectra # peaks');

plotX=25; plotY=15;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myPaths.reports,['Summary5_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r300');

% =========================================================================
% Report table
save(fullfile(myPaths.reports,['Summary_' myPaths.group '_' myPaths.visit '_' myPaths.task '_' myPaths.proctime]),'tableIssues');
writetable(tableIssues,fullfile(myPaths.reports,['Summary_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime '.xlsx']),"WriteMode","overwrite");

end