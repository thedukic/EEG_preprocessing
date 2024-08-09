function EEG = report_ICA(EEG)
%
% IC topoplots
%
% SDukic, August 2024

% Separate EEG
chanext = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EXT')).labels};
if ~isempty(chanext)
    EXT = pop_select(EEG,'channel',chanext);
    EEG = pop_select(EEG,'rmchannel',chanext);
end

% % Make sure IC activations are present
% % But also this is necessary if wICA was done!
% EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
% EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

% Evaluate only the first K = 20 ICs, they carry the most power and thus relevance
K = 20;
NICA = length(EEG.reject.gcompreject);
varICs = EEG.ALSUTRECHT.ica.varICs;

% Normalised var of all ICs
varAICsNorm = varICs./sum(varICs);
% Total var of the first K ICs
varKICs = sum(varAICsNorm(1:K));

% =========================================================================
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
    fprintf(EEG.ALSUTRECHT.subject.fid,'EEG data seems to be long enough for the ICA.\n');
end

% =========================================================================
% Plot the first 20 ICs
myCmap1 = brewermap(128,'*RdBu');

fh = figure;
th = tiledlayout(4,5);
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:20
    nexttile;
    topoplot(EEG.icawinv(:,i),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,i)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map');

    thisLabel = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)};
    if contains(thisLabel,'Brain')
        thisTitleColor = [0.1 0.8 0.2];
    elseif contains(thisLabel,'Other')
        thisTitleColor = [0 0 0];
    else
        thisTitleColor = [0.8 0.1 0.2];
    end

    % title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i),2))], [EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.ICLabel.cvec(i)} ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2))]});
    title({['ICA' num2str(i) ', Var = ' num2str(round(varAICsNorm(i)*100)) '%'], [thisLabel ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel.pvec(i),2))]},'Color',thisTitleColor);
    axis tight;
end

% Save
plotX=25; plotY=15;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ICs']),'-dtiff','-r300');
close(fh);

% =========================================================================
% % Plot IClabel artifact ICs
%
% fh = figure;
% th = tiledlayout('flow');
% th.TileSpacing = 'compact'; th.Padding = 'compact';
%
% for i = 1:length(EEG.ALSUTRECHT.ica.ICLabel.bics)
%     nexttile;
%     topoplot(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel.bics(i)),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel.bics(i))))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map');
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
%
% % =========================================================================
% % Plot extra artifact ICs
%
% % % Maximum number of bad ICs per each type
% % NLBL = length(EEG.ALSUTRECHT.ica.combi.clss);
% % maxBadIC = max(arrayfun(@(i) sum(EEG.ALSUTRECHT.ica.combi.cvec==i),1:NLBL));
% %
% % % Plot
% % fh = figure;
% % th = tiledlayout(NLBL,maxBadIC);
% % th.TileSpacing = 'compact'; th.Padding = 'compact';
% %
% % for i = 1:NLBL
% %     badICtypeIndx = find(EEG.ALSUTRECHT.ica.combi.cvec==i);
% %
% %     if ~isempty(badICtypeIndx)
% %         for j = 1:length(badICtypeIndx)
% %             nexttile(maxBadIC*(i-1)+j);
% %             thisIC = EEG.ALSUTRECHT.ica.combi.bics(badICtypeIndx(j));
% %             topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map');
% %             title({['ICA' num2str(thisIC)], [EEG.ALSUTRECHT.ica.combi.clss{i} ', R = ' num2str(round(EEG.ALSUTRECHT.ica.combi.corr(thisIC,i),2))]});
% %             axis tight;
% %         end
% %     end
% % end
%
% % Extra 1
% fh = figure;
% th = tiledlayout('flow');
% th.TileSpacing = 'compact'; th.Padding = 'compact';
%
% for i = 1:length(EEG.ALSUTRECHT.ica.extra1.bics)
%     nexttile;
%
%     thisIC = EEG.ALSUTRECHT.ica.extra1.bics(i);
%     thisClass = EEG.ALSUTRECHT.ica.extra1.cvec(i);
%     topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map');
%     title({['ICA' num2str(thisIC)], [EEG.ALSUTRECHT.ica.extra1.clss{thisClass} ', R = ' num2str(round(EEG.ALSUTRECHT.ica.extra1.corr(thisIC,thisClass),2))]});
%     axis tight;
% end
%
% % Save
% plotX=35; plotY=20;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_extra1_badICs']),'-dtiff','-r400');
% close(fh);
%
% % Extra 2
% fh = figure;
% th = tiledlayout('flow');
% th.TileSpacing = 'compact'; th.Padding = 'compact';
%
% for i = 1:length(EEG.ALSUTRECHT.ica.extra2.bics)
%     nexttile;
%     thisIC = EEG.ALSUTRECHT.ica.extra2.bics(i);
%     thisClass = EEG.ALSUTRECHT.ica.extra2.cvec(i);
%     topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map');
%     if thisClass == 1
%         measureLabel = 'Slope';
%         measureValue = EEG.ALSUTRECHT.ica.extra2.musleSlope(thisIC);
%     else
%         measureLabel = 'Pval';
%         measureValue = mean(EEG.ALSUTRECHT.ica.extra2.blinkPvals(thisIC,:));
%     end
%     title({['ICA' num2str(thisIC)], [EEG.ALSUTRECHT.ica.extra2.clss{thisClass} ', ' measureLabel ' = ' num2str(round(measureValue,2))]});
%     axis tight;
% end
%
% % Save
% plotX=35; plotY=20;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_extra2_badICs']),'-dtiff','-r400');
% close(fh);

% =========================================================================
% Plot all artifact ICs

fh = figure;
th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

myCmap2 = brewermap(4,'Set1');
for i = 1:length(EEG.ALSUTRECHT.ica.combi.bics)
    nexttile;
    thisIC = EEG.ALSUTRECHT.ica.combi.bics(i);
    topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map');
    title({['ICA' num2str(thisIC)], [EEG.ALSUTRECHT.ica.combi.lbls{i} ', ' EEG.ALSUTRECHT.ica.combi.meth{i} ' = ' num2str(round(EEG.ALSUTRECHT.ica.combi.prbs(i),2))]},'Color',myCmap2(EEG.ALSUTRECHT.ica.combi.method(i),:));
end

% Save
plotX=25; plotY=15;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_badICs']),'-dtiff','-r300');
close(fh);

% =========================================================================
% Plot wICA eye and muscle

ICsMostLikelyBlink  = EEG.ALSUTRECHT.ica.extra3.bics;
ICsMostLikelyMuscle = find(EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle);

% Eye
fh = figure;
th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:length(ICsMostLikelyBlink)
    nexttile;
    thisIC = ICsMostLikelyBlink(i);
    topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map');
    title(['ICA' num2str(thisIC),', Eye']);
end

% Save
plotX=20; plotY=10;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_badICsEye']),'-dtiff','-r300');
close(fh);

% Muscle
fh = figure;
th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:length(ICsMostLikelyMuscle)
    nexttile;
    thisIC = ICsMostLikelyMuscle(i);
    topoplot(EEG.icawinv(:,thisIC),EEG.chanlocs,'maplimits',max(abs(EEG.icawinv(:,thisIC)))*[-1 1],'headrad','rim','colormap',myCmap1,'whitebk','on','style','map');
    title(['ICA' num2str(thisIC),', Muscle']);
end

% Save
plotX=20; plotY=10;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_badICsMuscle']),'-dtiff','-r300');
close(fh);

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

fprintf(EEG.ALSUTRECHT.subject.fid,'Within the first %d ICs (power = %1.2f):\n', K,varKICs);
fprintf(EEG.ALSUTRECHT.subject.fid,'Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==1)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==2)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==3)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==4)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==5)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==6)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==7)*100));

fprintf('Within the first %d ICs (power = %1.2f):\n', K,varKICs);
fprintf('Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==1)*100));
fprintf('Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==2)*100));
fprintf('Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==3)*100));
fprintf('Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==4)*100));
fprintf('Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==5)*100));
fprintf('Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==6)*100));
fprintf('Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.combi.report(1:K)==7)*100));

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

% Merge
if ~isempty(chanext)
    EEG = merge_eeglabsets(EEG,EXT);
end

end
