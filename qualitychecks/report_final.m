function report_final(myfolders,rnum)
%
% TODO: trials rejected
%

fprintf('Generating a final report... This may take a while...\n');
if isnumeric(rnum), rnum = num2str(rnum); end

subjects = list_subjects(myfolders.preproc,[]);
subjects = subjects(contains(subjects,'ALS'));
NSUB = length(subjects);

chanlocs = readlocs('biosemi128_eeglab.ced');
chanlbls = {chanlocs.labels};
maskelec = zeros(length(chanlocs),NSUB);

% =========================================================================
N = NaN(NSUB,5);
Medianvoltageshift = NaN(length(chanlbls),NSUB);

for i = 1:NSUB
    load(fullfile(myfolders.preproc,subjects{i},[subjects{i} '_' myfolders.visit '_' myfolders.task '_cleandata_' rnum 'b.mat']),'preprocReport');

    maskelec(:,i) = double(ismember(chanlbls,preprocReport.badchaninfo.badElectrodes));
    N(i,1) = sum(maskelec(:,i));

    % Total possible
    % switch myfolders.task
    %     case {'RS','EO','EC'}
    %         % Assume 2s with 0.75 overlap
    %         N(i,2) = sum([preprocReport.eventinfo{:,3}])*4/2;
    %     otherwise
    %         % N(i,2) = sum([preprocReport.eventinfo{:,3}]);
    % end
    N(i,2) = preprocReport.issues_to_check.NumberTrials1;

    % Left after preproc
    N(i,3) = preprocReport.issues_to_check.NumberTrials2;
    % Left after trial rejection
    N(i,4) = preprocReport.issues_to_check.NumberTrials4;
    % Leftover EMG
    N(i,5) = preprocReport.issues_to_check.MuscleLeftovers;

    Medianvoltageshift(:,i) = preprocReport.issues_to_check.MedianvoltageshiftwithinepochFinal(1:128);
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
axis tight; title([myfolders.group ', N = ' num2str(NSUB)]);
hcb = colorbar;
hcb.Title.String = "%";

% Plot 2: Interpoalted channels per person
nexttile; hold on;
plot([0 NSUB],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,1)/128,'FaceColor',myCmap1(end/2,:)); ylim([0 1]); box on;

ylabel('Interp. chans (%)');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);

% Plot 3: EMG lefovers from preprocessing but might be less after trial rejection
nexttile; hold on;
plot([0 NSUB],[0.25 0.25],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,5),'FaceColor',myCmap3(12,:)); ylim([0 1]); box on;

xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);
ylabel('EMG leftover (%)');

% Plot 3
nexttile;
Ntmp = [N(:,4), N(:,3)-N(:,4), N(:,2)-N(:,3)];
bh = bar(Ntmp,'stacked');
bh(1).FaceColor = myCmap2(3,:);
bh(2).FaceColor = myCmap2(2,:);
bh(3).FaceColor = myCmap2(1,:);

ylim([0 max(sum(Ntmp,2))]);
ylabel('# trials');
xticks(1:NSUB); xticklabels(subjects); xlim([0 NSUB+1]);
% xlabel('Participant');

% Save
plotX=25; plotY=22;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myfolders.reports,['Summary1_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime]),'-dtiff','-r400');

% =========================================================================
% The following checks for participants who show outlying values for the median voltage shift within each epoch:
% The following detects outlier files in the median amount of their max-min
% voltage shift within an epoch, after adjusting for the fact that the data
% across all participants is likely to be positively skewed with a log transform.
MedianvoltageshiftwithinepochLogged=log10(Medianvoltageshift);
InterQuartileRange=iqr(MedianvoltageshiftwithinepochLogged,2);
Upper25 = prctile(MedianvoltageshiftwithinepochLogged,75,2);
Lower25 = prctile(MedianvoltageshiftwithinepochLogged,25,2);

% 75th% and 25th% +/- (1.5 x IQR) is the recommended outlier detection method,
% so this is used to recommend which participants to manually check
% However, I find this can be a bit too sensitive upon manual inspection,
% and that 1.75, 2, or even 2.5 can be a better threshold.
LowerBound=size(MedianvoltageshiftwithinepochLogged,1);
UpperBound=size(MedianvoltageshiftwithinepochLogged,1);
for i=1:size(MedianvoltageshiftwithinepochLogged,1)
    LowerBound(i,1)=Lower25(i,1)-(2*InterQuartileRange(i,1));
    UpperBound(i,1)=Upper25(i,1)+(2*InterQuartileRange(i,1));
end

VoltageShiftsTooLow=MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooLow=VoltageShiftsTooLow-LowerBound;
VoltageShiftsTooLow(VoltageShiftsTooLow>0)=0;
CumulativeSeverityOfAmplitudesBelowThreshold = sum(VoltageShiftsTooLow,1)';

VoltageShiftsTooHigh=MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooHigh=VoltageShiftsTooHigh-UpperBound;
VoltageShiftsTooHigh(VoltageShiftsTooHigh<0)=0;
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
print(fh,fullfile(myfolders.reports,['Summary2_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime]),'-dtiff','-r400');

% 3. Visualise
maskSubj = find(CumulativeSeverityOfAmplitudesAboveThreshold>0);

fh = figure;
th = tiledlayout(1,length(maskSubj));
th.TileSpacing = 'compact'; th.Padding = 'compact';
title(th,'Participants with high median voltage');

myClim = [mean(LowerBound), mean(UpperBound)];

for i = 1:length(maskSubj)
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
print(fh,fullfile(myfolders.reports,['Summary3_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime]),'-dtiff','-r400');

% OutlierParticipantsToManuallyCheck   = table(Participant_IDs', CumulativeSeverityOfAmplitudesBelowThreshold,CumulativeSeverityOfAmplitudesAboveThreshold);
% LoggedMedianVoltageShiftAcrossEpochs = array2table(MedianvoltageshiftwithinepochLogged);

end