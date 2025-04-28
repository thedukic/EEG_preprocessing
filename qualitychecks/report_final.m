function report_final(myPaths,subjects)
% =========================================================================
%
% Script for reporting on EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================

fprintf('\n==================================================================\n');
fprintf('Generating the final report for %s. This may take a bit...\n',myPaths.group);
fprintf('==================================================================\n');

% =========================================================================
NSUB = length(subjects);
chanlocs = readlocs('biosemi128_eeglab.ced');
chanlbls = {chanlocs.labels};

maskInterpElec1 = zeros(length(chanlocs),NSUB);
maskInterpElec2 = maskInterpElec1;

% Report folder
myPaths.reports = fullfile(myPaths.preproc,'reports');
if exist(myPaths.reports,'dir')~=7, mkdir(myPaths.reports); end

% =========================================================================
% Preprocessing stats
% =========================================================================
N = NaN(NSUB,6);
NCHN = length(chanlbls);
Medianvoltageshift = NaN(NCHN,NSUB);
CorrelationMatrices = NaN(NCHN+2,NCHN+2,NSUB);

fprintf('Loading %d datasets. Wait...\n',NSUB);
for i = 1:NSUB
    fileName1 = fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'a.mat']);
    fileName2 = fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat']);
    if exist(fileName2,"file") == 2
        % Load
        load(fileName2,'EEG');

        % Record warnings about potential issues
        EEG = report_issues(EEG);

        % Record warnings for all participants in a single table
        if i == 1
            tableIssues = struct2table(EEG.ALSUTRECHT.issues_to_check,'AsArray',true);
            tableIssues(2:NSUB,:) = cell2table(...
                [repmat({''},NSUB-1,1), repmat({NaN},NSUB-1,length(tableIssues.Properties.VariableNames)-1)], ...
                'VariableNames', tableIssues.Properties.VariableNames);
        else
            tableIssues(i,:) = struct2table(EEG.ALSUTRECHT.issues_to_check,'AsArray',true);
        end

        % Bad electrodes (whole interpolated)
        maskInterpElec1(:,i) = double(ismember(chanlbls,EEG.ALSUTRECHT.badchaninfo.badElectrodes));

        % Bad electrodes (trials interpolated)
        interpChanTrials = cat(1,EEG.ALSUTRECHT.epochRejections.InterpTrialInfo{:});
        interpChanTrials = unique(interpChanTrials);
        maskInterpElec2(interpChanTrials,i) = 1;

        % Number of interpolated channels
        N(i,1) = sum(maskInterpElec1(:,i)) / 128;
        % Number of interpolated trials
        N(i,2) = EEG.ALSUTRECHT.epochRejections.interpEpochs / length(EEG.ALSUTRECHT.epochRejections.InterpTrialInfo);

        % Epochs: Total possible
        % N(i,3) = sum([EEG.ALSUTRECHT.eventinfo{:,3}])*4/2; % RS
        N(i,3) = EEG.ALSUTRECHT.issues_to_check.NumberTrials1;
        % Epochs: Left after preproc1
        N(i,4) = EEG.ALSUTRECHT.issues_to_check.NumberTrials2;
        % Epochs: Left after preproc2
        N(i,5) = EEG.ALSUTRECHT.issues_to_check.NumberTrials3;

        % Leftover EMG
        % N(i,5) = EEG.ALSUTRECHT.leftovers.muscle1;  % after proc1
        N(i,6) = EEG.ALSUTRECHT.leftovers.muscle2;    % after proc2

        % Median voltage range
        Medianvoltageshift(:,i) = EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal(1:128);

        % Correlation matrix
        CorrelationMatrices(:,:,i) = EEG.ALSUTRECHT.chanCorr;

    elseif exist(fileName1,"file") == 2
        % Then, load the file from preproc1
        fprintf('%s does not have completely preprocessed %s data (only from part 1).\n', subjects{i}, myPaths.task);
        load(fileName1,'EEG');

        % Insert manually
        maskInterpElec1(:,i) = double(ismember(chanlbls,EEG.ALSUTRECHT.badchaninfo.badElectrodes));
        N(i,1) = sum(maskInterpElec1(:,i));
        N(i,2) = 0; % True
        N(i,3) = sum([EEG.ALSUTRECHT.eventinfo{:,3}]);
        N(i,4) = 0; % Not true but OK
        N(i,5) = 0; % True
        N(i,6) = 1; % Probably true

    else
        fprintf('%s does not have any preprocessed %s data.\n', subjects{i}, myPaths.task);
        N(i,:) = NaN;
    end
end

% Remove those with very bad or missing data
% -> they have 0 trials after preproc2
% -> they have NaNs
maskBad     = N(:,5) == 0;
maskMissing = isnan(N(:,1));
maskRemove  = maskBad | maskMissing;

Medianvoltageshift(:,maskRemove)    = [];
CorrelationMatrices(:,:,maskRemove) = [];
maskInterpElec1(:,maskRemove)       = [];
maskInterpElec2(:,maskRemove)       = [];

% % Double-check
% assert(~any(isnan(N(:))));

% ============================
% Plot
% ============================
fh = figure;
th = tiledlayout(5,2);
th.TileSpacing = 'compact'; th.Padding = 'compact';
title(th,[myPaths.group ', N = ' num2str(NSUB)]);

% Plot 1: Which channels are usually interpoalted
myCmap1 = brewermap(128,'RdPu');
myCmap2 = brewermap(3,'YlOrRd');
myCmap3 = brewermap(12,'Paired');

nexttile(1);
topoplot(mean(maskInterpElec1,2),chanlocs,'maplimits',[0 0.25],'headrad',0.5,'colormap',myCmap1,'whitebk','on','electrodes','on','style','map','shading','interp');
axis tight; title('Interp. chans: preproc1');
hcb = colorbar;
hcb.Title.String = "%";

nexttile(2);
topoplot(mean(maskInterpElec2,2),chanlocs,'maplimits',[0 0.7],'headrad',0.5,'colormap',myCmap1,'whitebk','on','electrodes','on','style','map','shading','interp');
axis tight; title('Interp. chans: preproc2');
hcb = colorbar;
hcb.Title.String = "%";

% Plot 2: Interpolated channels per person
nexttile([1 2]); hold on;
plot([0 NSUB+1],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,1),'FaceColor',myCmap1(end/2,:)); ylim([0 1]); box on;

ylabel('Interp. chans preproc1 (%)');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);

% Plot 2: Interpolated channels per person
nexttile([1 2]); hold on;
plot([0 NSUB+1],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,2),'FaceColor',myCmap1(end/2,:)); ylim([0 1]); box on;

ylabel('Interp. trials preproc2 (%)');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);

% Plot 3: EMG lefovers from preprocessing
nexttile([1 2]); hold on;
plot([0 NSUB+1],[0.25 0.25],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,6),'FaceColor',myCmap3(12,:)); ylim([0 1]); box on;

xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);
ylabel('EMG leftover (%)');

% Plot 4: Number of trials overview
nexttile([1 2]);
Ntmp = [N(:,5), N(:,4)-N(:,5), N(:,3)-N(:,4)];
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
% Voltage range/shifts across channels
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

% ============================
% Plot
% ============================
% Plot 1
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

% Plot 2
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
        maskChan = VoltageShiftsTooHigh(:,maskSubj(i)) > 0;
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
% Channel correlation
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
% Power spectra deviations in gamma-band
% -> using stat. measusre to detect large spread across channels
% =========================================================================
warning('The cutoff for spectra spread is tested only for resting-state. It might not be appropriate for other tasks.');
cutoffKurtosis = Inf; % dont use
cutoffCoefVar  = 1;   % RS
cutoffPeak     = 5;

cnt = 0;
psdSpreadsAll1 = NaN(NSUB,1);
psdSpreadsAll2 = NaN(NSUB,1);
psdPeaksAll    = NaN(NSUB,1);
psdPromsAll    = NaN(NSUB,1);
psdSpreadAll   = NaN(NSUB,NCHN);

fprintf('Loading %d datasets... Wait...\n',NSUB);
for i = 1:NSUB
    % Do this only if data is fully preprocessed
    fileName = fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat']);
    if exist(fileName,"file") == 2
        % Load
        load(fileName,'EEG');

        % Estimate the spectra
        [psdspectra, freq, chaneeg] = estimate_power(EEG,'freport');
        % psdspectra = log10(psdspectra);
        % psdspectraNorm = psdspectra ./ sum(psdspectra(freq>0,:),1);

        % ============================
        % 1. Calculate measure of power spread across EEG electrodes
        % ============================
        maskFreq  = freq>25 & freq<45;
        % psdSpread = psdspectra ./ sum(psdspectra,1);
        psdSpread = mean(psdspectra(maskFreq,:),1);

        % Spread
        psdSpreadsAll1(i) = kurtosis(psdSpread);
        psdSpreadsAll2(i) = std(psdSpread) ./ mean(psdSpread);

        % Use for clustering/outlier detection (later)
        psdSpreadAll(i,:) = psdSpread;

        % ============================
        % 2. Find the number of peaks in the average EEG spectrum
        % ============================
        [locs, psdPromsAll(i)] = find_peaks(psdspectra,freq);
        psdPeaksAll(i) = length(locs);

        % ============================
        % Plot and output if it exceeds the cut-off
        % ============================
        if psdSpreadsAll1(i) > cutoffKurtosis || psdSpreadsAll2(i) > cutoffCoefVar || psdPeaksAll(i) > cutoffPeak
            fprintf('%d. %s : Kurtosis = %1.3f, CV = %1.3f, Peaks = %d\n', i, subjects{i}, psdSpreadsAll1(i), psdSpreadsAll2(i), psdPeaksAll(i));

            cnt = cnt+1;
            if cnt == 1
                dataCmap = brewermap(sum(chaneeg),'BrBG');

                fh = figure;
                th = tiledlayout("flow");
                th.TileSpacing = 'compact'; th.Padding = 'compact';
            end

            th = nexttile;
            hold on; box off;

            % psdspectra  = log10(psdspectra);
            % dataClim    = 1.1*[min(psdspectra(:)), max(psdspectra(:))];
            % dataClim(1) = floor(dataClim(1));
            % dataClim(2) = ceil(dataClim(2));
            plot(freq,log10(psdspectra),'LineWidth',1.1);
            colororder(th,dataCmap);

            for j = 1:psdPeaksAll(i)
                plot(locs(j) * ones(1,2),[-4 3],'LineWidth',1,'Color',0.2*ones(1,3));
                % scatter(locs,psdspectra(any(freq==locs',2)),'filled','*','MarkerFaceColor',0.2*ones(1,3));
            end

            xticks(unique([locs' 0 5 10 15 20 25 30 50]));
            xlim([freq(1), 60]); ylim([-4 3]); % ylim(dataClim);
            % title({[subjects{i},', group: ',myPaths.group, ', visit: ', myPaths.visit, ', task: ', myPaths.task], ['Spread: ', num2str(psdSpreadsAll1(i), '%.3f') ', Peaks: ', num2str(psdPeaksAll(i), '%d') ]});
            title({subjects{i}, ['Kurtosis: ', num2str(psdSpreadsAll1(i), '%1.1f') ', CV: ', num2str(psdSpreadsAll2(i), '%.2f') ', Peaks: ', num2str(psdPeaksAll(i), '%d') ]});

            pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('log_{10}(Power) (a.u.)'); drawnow;
        end
    end
end

% Save
if cnt > 0
    plotX=25; plotY=25;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(myPaths.reports,['Summary4_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r300');
    % close(fh);
end

% ============================
% Plot some additional info
% ============================
% Plot the distributions of the power spectra characteristics
% -> Power spectra spread, power spectra peaks
fh = figure;
th = tiledlayout(2,2);
th.TileSpacing = 'compact'; th.Padding = 'compact';
nexttile;
hb = histogram(psdSpreadsAll1,10);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power spectra: kurtosis');
myXlim = xlim; xlim([0 myXlim(2)]);
nexttile;
hb = histogram(psdSpreadsAll2,10);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power spectra: coefficient of variation');
myXlim = xlim; xlim([0 myXlim(2)]);
nexttile;
hb = histogram(psdPeaksAll);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power peaks');
myXlim = xlim; xlim([0 myXlim(2)]);
nexttile;
hb = histogram(psdPromsAll,10);
hb.FaceColor = 0.7*ones(1,3);
pbaspect([1.618 1 1]); xlabel('Power peak prominance');
myXlim = xlim; xlim([0 myXlim(2)]);
% nexttile; plot(psdPromsAll);

% Save
plotX=20; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myPaths.reports,['Summary5_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime]),'-dtiff','-r300');

% =========================================================================
% Power spectra deviations in gamma-band
% -> using DBSCAN clustering and Mahalanobis distance
% =========================================================================
% Normalise
psdSpreadAll = psdSpreadAll';
psdSpreadAll = psdSpreadAll - min(psdSpreadAll(:));
psdSpreadAll = psdSpreadAll / max(psdSpreadAll(:));
psdSpreadAll = bsxfun(@minus, psdSpreadAll, mean(psdSpreadAll, 1));
psdSpreadAll = psdSpreadAll';

% figure; histogram(psdSpreadAll(:));

% =========================================
% PCA
[coeff, score, latent, tsquared, explained, mu] = pca(psdSpreadAll);

explained2 = explained ./ sum(explained);
Nkeep = 3;

fh = figure;
th = tiledlayout(1,3);
th.TileSpacing = 'compact'; th.Padding = 'compact';
nexttile; bar(explained2); title(['PCA1-' num2str(Nkeep) ' = ' num2str(round(sum(explained2(1:Nkeep)),2))]); pbaspect([1.618 1 1]);
nexttile;
scatter(score(:,1),score(:,2),'filled');
xlabel('PCA1'); ylabel('PCA2'); pbaspect([1.618 1 1]);
nexttile;
scatter3(score(:,1),score(:,2),score(:,3),'filled');
xlabel('PCA1'); ylabel('PCA2'); zlabel('PCA3'); pbaspect([1.618 1 1]);

% =========================================
% DBSCAN
X = score(:,1:Nkeep);
clustRes = dbscan(X,0.25,3);

Nclust = length(unique(clustRes));
myCmap = brewermap(Nclust,'Set2');

assert(Nkeep == 3);
% assert(Nclust == 2);

figure;
th = tiledlayout(1,3);
th.TileSpacing = 'compact'; th.Padding = 'compact';

nexttile; hold on;
% gscatter(score(:,1),score(:,2),clustRes);
% xlabel('PCA1'); ylabel('PCA2'); pbaspect([1.618 1 1]);

clustRes2 = clustRes;
clustRes2(clustRes2>0) = clustRes2(clustRes2>0)+1;
clustRes2(clustRes2==-1) = 1;

for i = 1:Nclust
    MyPlotData = X(clustRes2==i, 1:3);
    plot3(MyPlotData(:,1),MyPlotData(:,2),MyPlotData(:,3),'Color',myCmap(i,:),'LineStyle','none','Marker','.','MarkerSize',10);
    ThisBoundary = boundary(MyPlotData);
    if size(MyPlotData,1)==2
        plot3(MyPlotData(:,1),MyPlotData(:,2),MyPlotData(:,3),'Color',myCmap(i,:),'LineWidth',2);
    elseif size(MyPlotData,1)==3
        patch(MyPlotData(:,1),MyPlotData(:,2),MyPlotData(:,3),myCmap(i,:),'Facecolor',myCmap(i,:),'Edgecolor',myCmap(i,:),'FaceAlpha',0.6)
    else
        trisurf(ThisBoundary,MyPlotData(:,1),MyPlotData(:,2),MyPlotData(:,3),'Facecolor',myCmap(i,:),'Edgecolor',myCmap(i,:),'FaceAlpha',0.6)
    end
end
grid on; box on; axis vis3d; view(-30,20);
xlabel('PCA1'); ylabel('PCA2'); zlabel('PCA3');

m = clustRes==-1;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],['Outliers: ' num2str(sum(m))],nexttile); colorbar;
m = clustRes>-1;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],['Rest: ' num2str(sum(m))],nexttile); colorbar;

% =========================================
% Mahalanobis
X = score(:,1:Nkeep);

% Calculate the Mahalanobis distance for each subject
D_M = mahal(X, X);

% Identify potential outliers
threshold = chi2inv(0.975, size(X, 2)); % 97.5% confidence interval
outliers = find(D_M > threshold);

disp('Potential outliers:');
disp(subjects(outliers)');

fh = figure;
th = tiledlayout(2, ceil(length(outliers)/2));
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:length(outliers)
    tmp = psdSpreadAll(outliers(i),:);
    mytopoplot(tmp,[],subjects{outliers(i)},nexttile); colorbar;
end

title(th,{'Participants with high Mahalanobis distance'});

% =========================================================================
% Report table
save(fullfile(myPaths.reports,['Summary_' myPaths.group '_' myPaths.visit '_' myPaths.task '_' myPaths.proctime]),'tableIssues');
writetable(tableIssues,fullfile(myPaths.reports,['Summary_' myPaths.group '_' myPaths.visit '_' myPaths.task  '_' myPaths.proctime '.xlsx']),"WriteMode","overwrite");

end