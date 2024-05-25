function report_final(myfolders,rnum)
%
% TODO: trials rejected
%

fprintf('Generating a final report... This may take a while...\n');
if isnumeric(rnum), rnum = num2str(rnum); end

subjects = list_subjects(myfolders.preproc,[]);
NSUB = length(subjects);

chanlocs = readlocs('biosemi128_eeglab.ced');
chanlbls = {chanlocs.labels};
maskelec = zeros(length(chanlocs),NSUB);

N = NaN(NSUB,4);
for i = 1:NSUB
    load(fullfile(myfolders.preproc,subjects{i},[subjects{i} '_' myfolders.visit '_' myfolders.task '_cleandata_' rnum '.mat']),'EEG');

    maskelec(:,i) = double(ismember(chanlbls,EEG.ALSUTRECHT.badchaninfo.badElectrodes));
    N(i,1) = sum(maskelec(:,i));

    % Total possible
    switch myfolders.task
        case {'RS','EO','EC'}
            % Assume 2s with 0.75 overlap
            N(i,2) = sum([EEG.ALSUTRECHT.eventinfo{:,3}])*4/2;
        otherwise
            N(i,2) = sum([EEG.ALSUTRECHT.eventinfo{:,3}]);
    end
    % Left after preproc
    N(i,3) = size(EEG.data,3);
    % Leftover EMG
    N(i,4) = EEG.ALSUTRECHT.issues_to_check.MuscleLeftovers;
end

% Plot
fh = figure;
th = tiledlayout(4,1);
th.TileSpacing = 'compact'; th.Padding = 'compact';

% Plot 1
myCmap = brewermap(256,'RdPu');

nexttile(1);
topoplot(mean(maskelec,2),chanlocs,'maplimits',[0 0.25],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','on','style','map');
axis tight; title(myfolders.group);
hcb = colorbar;
hcb.Title.String = "%";

% Plot 2
nexttile; hold on;
plot([0 NSUB],[0.15 0.15],'LineWidth',1.2,'Color',0.6*ones(1,3));
bar(N(:,1)/128); ylim([0 0.25]); box on;
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
print(fh,fullfile(myfolders.preproc,['Summary_' myfolders.group '_' myfolders.visit '_' myfolders.task  '_' myfolders.proctime]),'-dtiff','-r400');

end