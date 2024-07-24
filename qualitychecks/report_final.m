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
N = NaN(NSUB,4);
Medianvoltageshift = NaN(length(chanlbls),NSUB);

for i = 1:NSUB
    load(fullfile(myfolders.preproc,subjects{i},[subjects{i} '_' myfolders.visit '_' myfolders.task '_cleandata_' rnum '.mat']),'preprocReport');

    maskelec(:,i) = double(ismember(chanlbls,preprocReport.badchaninfo.badElectrodes));
    N(i,1) = sum(maskelec(:,i));

    % Total possible
    switch myfolders.task
        case {'RS','EO','EC'}
            % Assume 2s with 0.75 overlap
            N(i,2) = sum([preprocReport.eventinfo{:,3}])*4/2;
        otherwise
            N(i,2) = sum([preprocReport.eventinfo{:,3}]);
    end
    % Left after preproc
    N(i,3) = preprocReport.issues_to_check.NumberTrials2;
    % Leftover EMG
    N(i,4) = preprocReport.issues_to_check.MuscleLeftovers;

    Medianvoltageshift(:,i) = preprocReport.issues_to_check.Medianvoltageshiftwithinepoch(1:128);
end

% =========================================================================
% Plot
fh = figure;
th = tiledlayout(4,1);
th.TileSpacing = 'compact'; th.Padding = 'compact';

% Plot 1
myCmap = brewermap(256,'RdPu');

nexttile(1);
topoplot(mean(maskelec,2),chanlocs,'maplimits',[0 0.25],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map');
axis tight; title([myfolders.group ', N = ' num2str(NSUB)]);
hcb = colorbar;
hcb.Title.String = "%";

% Plot 2
nexttile; hold on;
plot([0 NSUB],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,1)/128); ylim([0 0.35]); box on;
ylabel('Interpolated channels (%)');

% Plot 3
nexttile;
Ntmp = [N(:,3), N(:,2)-N(:,3)];
bar(Ntmp,'stacked');
ylabel('Number of trials');

% Plot 4
nexttile; hold on;
plot([0 NSUB],[0.25 0.25],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,4)); ylim([0 1]); box on;
ylabel('EMG leftover');
xlabel('Participant');

% Save
plotX=18; plotY=14;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(myfolders.reports,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime]),'-dtiff','-r400');

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
VoltageShiftsTooLow(0<VoltageShiftsTooLow)=0;
CumulativeSeverityOfAmplitudesBelowThreshold = sum(VoltageShiftsTooLow,1)';

VoltageShiftsTooHigh=MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooHigh=VoltageShiftsTooHigh-UpperBound;
VoltageShiftsTooHigh(0>VoltageShiftsTooHigh)=0;
CumulativeSeverityOfAmplitudesAboveThreshold = sum(VoltageShiftsTooHigh,1)';

% Plot:
fh = figure; hold on;
plot(LowerBound,'LineWidth',1.5,'Color','k'); 
plot(UpperBound,'LineWidth',1.5,'Color','k');
plot(MedianvoltageshiftwithinepochLogged,'LineWidth',1.2);

xticks(1:128); xticklabels({chanlocs.labels}); xlim([1 128]);
legend('LowerBound', 'UpperBound');

set(gca,'ColorOrder',[0 0 0; 0 0 0; brewermap(NSUB,'BuGn')]);

% OutlierParticipantsToManuallyCheck   = table(Participant_IDs', CumulativeSeverityOfAmplitudesBelowThreshold,CumulativeSeverityOfAmplitudesAboveThreshold);
% LoggedMedianVoltageShiftAcrossEpochs = array2table(MedianvoltageshiftwithinepochLogged);

end