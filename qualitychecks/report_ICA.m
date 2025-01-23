function EEG = report_ICA(EEG)
%
% IC topoplots
%

% =========================================================================
% Normalised variance of all ICs
varICs = EEG.ALSUTRECHT.ica.varICs;
varAICsNorm = varICs./sum(varICs);

% Evaluate only the first K ICs,
% They carry the most power and thus relevance
NICA = length(EEG.reject.gcompreject);
NICArel = 20;

% Total variance of the first K ICs
varKICs = sum(varAICsNorm(1:NICArel));

% =========================================================================
% Check if the data was enough for ICA

% EEG available:
% MMN  3*7     ~ 21 min
% SART 3*5     ~ 15 min
% RS   2x3x2   ~ 12 min
% MT   7+3+7   ~ 17 min
%
% Minimum EEG needed:
% 128 elecs ^2*30 /256/60 ~ 32  min
% 50 PCs    ^2*30 /256/60 ~ 5.0 min
% 70 PCs    ^2*30 /256/60 ~ 9.5 min
% See: https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline#What_is_the_minimum_length_of_data_to_perform_ICA.3F_.2807.2F04.2F2022_added.29

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'ICA\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
if 30*(NICA^2)>prod(size(EEG.data,[2 3]))
    EEG.ALSUTRECHT.ica.DataLengthForValidICA = 'NOK';
    warning('EEG data might be too short for ICA!');
    fprintf(EEG.ALSUTRECHT.subject.fid,'EEG data might be too short for the ICA!\n');
else
    EEG.ALSUTRECHT.ica.DataLengthForValidICA = 'OK';
    fprintf('\nEEG data is long enough for the ICA.\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'EEG data is long enough for the ICA.\n');
end

% =========================================================================
% Plot the first 20 ICs + bad ICs
myCmap1 = brewermap(128,'*RdBu');
myCmap2 = brewermap(4,'Set1');

% Extract bad ICs for plotting
% ICsforRemoval = [];
ICsforRemoval        = find(EEG.ALSUTRECHT.ica.ICsforRemoval);
ICsMostLikelyMuscle  = find(EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle2);
ICsMostLikelyChannel = find(EEG.ALSUTRECHT.ica.ICsMostLikelyChannel2);
ICsMostLikelyComplex = find(EEG.ALSUTRECHT.ica.ICsMostLikelyComplex);

NICAtmp = length([ICsforRemoval; ICsMostLikelyMuscle; ICsMostLikelyChannel; ICsMostLikelyComplex]);

% How many rows are needed for bad ICs
% We plot the first 24 ICs
NCOL = 8;
NROW1 = 3;
NICAgood = NCOL*NROW1;

NROW2 = ceil(NICAtmp/NCOL);

% Total rows plus 1 for the gap
NROW = NROW1 + NROW2 + 1;

fh = figure;
th = tiledlayout(NROW,NCOL);
th.TileSpacing = 'tight'; th.Padding = 'tight';

% 1. Plot the first N ICs
for i = 1:NICAgood
    nexttile;
    topoplot(EEG.icawinv(:,i),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,i)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','interp');

    thisLabel = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)};
    thisPval  = round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2);
    if contains(thisLabel,'Brain')
        thisTitleColor = [0.1 0.8 0.2];
    elseif contains(thisLabel,'Other')
        thisTitleColor = [0 0 0];
    else
        thisTitleColor = [0.8 0.1 0.2];
    end

    % title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i),2))], [EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)} ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2))]});
    title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i)*100)) '%'], [thisLabel ', P = ' num2str(thisPval)]},'Color',thisTitleColor);
    axis tight;
end

% 2. Plot bad ICs
NSTART = NICAgood + NCOL;
cnt = 0;
for i = 1:4
    switch i
        case 1
            theseICs  = ICsforRemoval;
            thisLabel = 'Removed';
        case 2
            theseICs  = ICsMostLikelyMuscle;
            thisLabel = 'Muscle';
        case 3
            theseICs  = ICsMostLikelyChannel;
            thisLabel = 'Channel';
        case 4
            theseICs = ICsMostLikelyComplex;
            thisLabel = 'Complex';
    end
    for j = 1:length(theseICs)
        cnt = cnt+1;
        nexttile(NSTART+cnt);

        thisIC = theseICs(j);
        switch i
            case 1
                % thisLabelTmp = EEG.ALSUTRECHT.ica.combi.lbls(EEG.ALSUTRECHT.ica.combi.bics==thisIC);
                % if length(thisLabelTmp)>1, thisLabelTmp = {strjoin(thisLabelTmp,'/')}; end
                thisLabelTmp = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.combi.report(thisIC)};
            case {2,3,4}
                thisLabelTmp = thisLabel;
        end

        topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','interp');
        title(['ICA' num2str(thisIC) ', ' thisLabelTmp],'Color',myCmap2(i,:));
    end
end

% Save
plotX=40; plotY=25;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_overviewICs']),'-dtiff','-r300');
close(fh);

% =========================================================================
% Plot all artifact ICs

% fh = figure;
% th = tiledlayout('flow');
% th.TileSpacing = 'tight'; th.Padding = 'tight';
%
% for i = 1:length(EEG.ALSUTRECHT.ica.combi.bics)
%     nexttile;
%     thisIC = EEG.ALSUTRECHT.ica.combi.bics(i);
%     topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','interp');
%     title({['ICA' num2str(thisIC)], [EEG.ALSUTRECHT.ica.combi.lbls{i} ', ' EEG.ALSUTRECHT.ica.combi.meth{i} ' = ' num2str(round(EEG.ALSUTRECHT.ica.combi.prbs(i),2))]},'Color',myCmap2(EEG.ALSUTRECHT.ica.combi.method(i),:));
% end
%
% % Save
% plotX=25; plotY=15;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_allbadICs']),'-dtiff','-r300');
% close(fh);

% =========================================================================
% % Plot the first 20 ICs
% myCmap1 = brewermap(128,'*RdBu');
%
% fh = figure;
% th = tiledlayout(4,5);
% th.TileSpacing = 'tight'; th.Padding = 'tight';
%
% for i = 1:20
%     nexttile;
%     topoplot(EEG.icawinv(:,i),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,i)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','interp');
%
%     thisLabel = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)};
%     thisPval  = round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2);
%     if contains(thisLabel,'Brain')
%         thisTitleColor = [0.1 0.8 0.2];
%     elseif contains(thisLabel,'Other')
%         thisTitleColor = [0 0 0];
%     else
%         thisTitleColor = [0.8 0.1 0.2];
%     end
%
%     % title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i),2))], [EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)} ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2))]});
%     title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i)*100)) '%'], [thisLabel ', P = ' num2str(thisPval)]},'Color',thisTitleColor);
%     axis tight;
% end
%
% % Save
% plotX=25; plotY=15;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ICs']),'-dtiff','-r300');
% close(fh);

% =========================================================================
% % Plot IClabel artifact ICs
%
% fh = figure;
% th = tiledlayout('flow');
% th.TileSpacing = 'tight'; th.Padding = 'tight';
%
% for i = 1:length(EEG.ALSUTRECHT.ica.ICLabel.bics)
%     nexttile;
%     topoplot(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel.bics(i)),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel.bics(i))))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map','shading','interp');
%     title({['ICA' num2str(EEG.ALSUTRECHT.ica.ICLabel.bics(i))], [EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics(i))} ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel.pvec(EEG.ALSUTRECHT.ica.ICLabel.bics(i)),2))]});
%     axis tight;
% end
%
% % Save
% plotX=35; plotY=20;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ICLabel_badICs']),'-dtiff','-r400');
% close(fh);

% =========================================================================
% % Plot wICs
% ICsforwICA           = EEG.ALSUTRECHT.ica.ICsforwICA;
% ICsMostLikelyBlink   = EEG.ALSUTRECHT.ica.ICsMostLikelyBlink;
% ICsMostLikelyMuscle  = EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle;
% ICsMostLikelyComplex = EEG.ALSUTRECHT.ica.ICsMostLikelyComplex;
%
% NICAtmp = sum(ICsforwICA|ICsMostLikelyBlink|ICsMostLikelyMuscle|ICsMostLikelyComplex);
%
% ICsforwICA           = find(ICsforwICA);
% ICsMostLikelyMuscle  = find(ICsMostLikelyMuscle);
% ICsMostLikelyComplex = find(ICsMostLikelyComplex);
%
% fh = figure;
% if NICAtmp<=5
%     th = tiledlayout(1,5);
% else
%     th = tiledlayout('flow');
% end
% th.TileSpacing = 'tight'; th.Padding = 'tight';
%
% for i = 1:3
%     switch i
%         case 1
%             theseICs  = ICsforwICA;
%             thisLabel = 'wICA';
%         case 2
%             theseICs  = ICsMostLikelyMuscle;
%             thisLabel = 'Muscle';
%         case 3
%             theseICs = ICsMostLikelyComplex;
%             thisLabel = 'Complex';
%     end
%     for j = 1:length(theseICs)
%         nexttile;
%         thisIC = theseICs(j);
%
%         switch i
%             case 1
%                 thisLabelTmp = EEG.ALSUTRECHT.ica.combi.lbls(EEG.ALSUTRECHT.ica.combi.bics==thisIC);
%                 if length(thisLabelTmp)>1, thisLabelTmp = {strjoin(thisLabelTmp,'/')}; end
%             case {2,3}
%                 thisLabelTmp = {thisLabel};
%         end
%
%         topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map','shading','interp');
%         title(['ICA' num2str(thisIC) ', ' thisLabelTmp{1}],'Color',myCmap2(i,:));
%     end
% end
%
% % Save
% plotX=25; plotY=15;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_badwICs']),'-dtiff','-r300');
% close(fh);

% =========================================================================
% Estimates of "good/successful" ICA:
%
% https://sccn.ucsd.edu/pipermail/eeglablist/2020/015096.html
% https://doi.org/10.1016/j.eplepsyres.2021.106809
% In Figure 4,
% During the awake state, Brain 53%, Muscle 12%, Eye 9%, Channel Noise <1%, Line noise <1%, Heart < 1%, Other 24%
% During the sleep state, Brain 74%, Other 23%, everything else < 2%
% This is for the case of 19 channels.
% Though this cross-number-of-channel test is not official,
% my impression is that the EEG data with standard quality seem to
% show 50-55 of Brain class rate regardless of the number of channels.
% (see also: eeglablist Digest, Vol 220, Issue 21)

fprintf(EEG.ALSUTRECHT.subject.fid,'Within the first %d ICs (power = %1.2f):\n', NICArel,varKICs);
fprintf(EEG.ALSUTRECHT.subject.fid,'Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==1)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==2)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==3)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==4)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==5)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==6)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==7)*100));

fprintf('Within the first %d ICs (power = %1.2f):\n', NICArel,varKICs);
fprintf('Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==1)*100));
fprintf('Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==2)*100));
fprintf('Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==3)*100));
fprintf('Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==4)*100));
fprintf('Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==5)*100));
fprintf('Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==6)*100));
fprintf('Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:NICArel)==7)*100));

% RELAX reporting:
I = EEG.ALSUTRECHT.ica.combi.report;
ICsMostLikelyBrain        = (I==1)';
ICsMostLikelyMuscle       = (I==2)';
ICsMostLikelyEye          = (I==3)';
ICsMostLikelyHeart        = (I==4)';
ICsMostLikelyLineNoise    = (I==5)';
ICsMostLikelyChannelNoise = (I==6)';
ICsMostLikelyOther        = (I==7)';

BrainVariance    = sum(abs(varICs(ICsMostLikelyBrain)));
ArtifactVariance = sum(abs(varICs(~ICsMostLikelyBrain)));
TotalVariance    = BrainVariance+ArtifactVariance;

MuscleVariance       = sum(abs(varICs(ICsMostLikelyMuscle)));
EyeVariance          = sum(abs(varICs(ICsMostLikelyEye)));
HeartVariance        = sum(abs(varICs(ICsMostLikelyHeart)));
LineNoiseVariance    = sum(abs(varICs(ICsMostLikelyLineNoise)));
ChannelNoiseVariance = sum(abs(varICs(ICsMostLikelyChannelNoise)));
OtherVariance        = sum(abs(varICs(ICsMostLikelyOther)));

EEG.ALSUTRECHT.ica.ProportionVariance_was_BrainICs        = BrainVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_MuscleICs       = MuscleVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_EyeICs          = EyeVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_HeartICs        = HeartVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_LineNoiseICs    = LineNoiseVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_ChannelNoiseICs = ChannelNoiseVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_OtherICs        = OtherVariance/TotalVariance;

% Remove (not needed)
EEG.icaact = [];

end