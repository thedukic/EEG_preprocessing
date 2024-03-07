function EEG = report_ica(EEG)
%
% IC topoplots
%
% SDukic, Feb 2024

% Make sure IC activations are present
if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

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
NICA = length(EEG.reject.gcompreject);

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
% Plot IClabel artifact ICs
% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel_bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel_clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel_pvec, EEG.ALSUTRECHT.ica.ICLabel_cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

myCmap = brewermap(128,'*RdBu');

fh = figure;
th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:length(EEG.ALSUTRECHT.ica.ICLabel_bics)
    nexttile;
    topoplot(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel_bics(i)),EEG.chanlocs,'maplimits',0.95*max(abs(EEG.icawinv(:,EEG.ALSUTRECHT.ica.ICLabel_bics(i))))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map');
    title({['ICA' num2str(EEG.ALSUTRECHT.ica.ICLabel_bics(i))], [EEG.ALSUTRECHT.ica.ICLabel_clss{EEG.ALSUTRECHT.ica.ICLabel_cvec(EEG.ALSUTRECHT.ica.ICLabel_bics(i))} ', P = ' num2str(round(EEG.ALSUTRECHT.ica.ICLabel_pvec(EEG.ALSUTRECHT.ica.ICLabel_bics(i)),2))]});
end

% Save
plotX=35; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ICLabel_badICs']),'-dtiff','-r400');
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
%

K = 20; % Evaluate only the first 20 ICs s they carry the most power and thus relevance
fprintf(EEG.ALSUTRECHT.subject.fid,'Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==1)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==2)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==3)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==4)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==5)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==6)*100));
fprintf(EEG.ALSUTRECHT.subject.fid,'Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==7)*100));

fprintf('Brain   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==1)*100));
fprintf('Muscle  components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==2)*100));
fprintf('Eye     components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==3)*100));
fprintf('Heart   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==4)*100));
fprintf('Line    components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==5)*100));
fprintf('Channel components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==6)*100));
fprintf('Other   components: %2.0f%%\n', round(mean(EEG.ALSUTRECHT.ica.ICLabel_cvec(1:K)==7)*100));

% RELAX reporting:
I = EEG.ALSUTRECHT.ica.ICLabel_cvec;
ICsMostLikelyBrain        = (I==1)';
ICsMostLikelyMuscle       = (I==2)';
ICsMostLikelyEye          = (I==3)';
ICsMostLikelyHeart        = (I==4)';
ICsMostLikelyLineNoise    = (I==5)';
ICsMostLikelyChannelNoise = (I==6)';
ICsMostLikelyOther        = (I==7)';

varianceWav = NaN(NICA,1);
for i = 1:NICA
    [~, varianceWav(i)] = compvar(EEG.data,EEG.icaact,EEG.icawinv,i);
end

fprintf(EEG.ALSUTRECHT.subject.fid,'Total power of the first %d ICs: %2.0f\n', K,sum(varianceWav(1:K)));
fprintf('Total power of the first %d ICs: %2.0f\n', K,sum(varianceWav(1:K)));

BrainVariance    = sum(abs(varianceWav(ICsMostLikelyBrain)));
ArtifactVariance = sum(abs(varianceWav(~ICsMostLikelyBrain)));
TotalVariance    = BrainVariance+ArtifactVariance;

MuscleVariance       = sum(abs(varianceWav(ICsMostLikelyMuscle)));
EyeVariance          = sum(abs(varianceWav(ICsMostLikelyEye)));
HeartVariance        = sum(abs(varianceWav(ICsMostLikelyHeart)));
LineNoiseVariance    = sum(abs(varianceWav(ICsMostLikelyLineNoise)));
ChannelNoiseVariance = sum(abs(varianceWav(ICsMostLikelyChannelNoise)));
OtherVariance        = sum(abs(varianceWav(ICsMostLikelyOther)));

EEG.ALSUTRECHT.ica.ProportionVariance_was_BrainICs        = BrainVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_MuscleICs       = MuscleVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_EyeICs          = EyeVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_HeartICs        = HeartVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_LineNoiseICs    = LineNoiseVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_ChannelNoiseICs = ChannelNoiseVariance/TotalVariance;
EEG.ALSUTRECHT.ica.ProportionVariance_was_OtherICs        = OtherVariance/TotalVariance;

end