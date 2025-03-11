function EEG = remove_badcomponents(EEG,cfg)
%
% https://www.biorxiv.org/content/10.1101/2024.06.06.597688v1.full.pdf
% SDukic, January 2025
% =========================================================================

fprintf('\n================================\n');
fprintf('Dealing with bad ICs\n');
fprintf('================================\n');

% Muscle ICs will be filtered such that they have only freq above this one
muscleFreq = cfg.ica.emgfilt; % [Hz]

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
dataICs = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
dataICs = reshape(dataICs, size(dataICs,1), EEG.pnts, EEG.trials);
dataICs = double(dataICs);

[NICA, NTPT] = size(dataICs);
assert(NTPT == EEG.pnts);
assert(ismatrix(dataICs));

% Allocate
artifactComponents = zeros(NICA,NTPT);

% Extract bad ICs
ICsMostLikelyEye     = EEG.ALSUTRECHT.ica.ICsMostLikelyEye;
ICsMostLikelyHeart   = EEG.ALSUTRECHT.ica.ICsMostLikelyHeart;
ICsMostLikelyComplex = EEG.ALSUTRECHT.ica.ICsMostLikelyComplex;
ICsMostLikelyMuscle  = EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle;
ICsMostLikelyChannel = EEG.ALSUTRECHT.ica.ICsMostLikelyChannel;

% Maybe we dont care about ICs with low variance
Nfix = 30;
fprintf('\nMuscle ICs above the %dth IC are not corrected.\n',Nfix);
fprintf('Channel ICs above the %dth IC are not corrected.\n\n',Nfix);
ICsMostLikelyMuscle((Nfix+1):end)  = false;
ICsMostLikelyChannel((Nfix+1):end) = false;

% *Remove completely ICs
ICsforRemoval = false(NICA,1);
ICsforRemoval(ICsMostLikelyEye | ICsMostLikelyHeart  | ICsMostLikelyComplex) = true;

% Report
fprintf('Eye ICs,     N = %d (removed completely).\n',sum(ICsMostLikelyEye));
fprintf('Heart ICs,   N = %d (removed completely).\n',sum(ICsMostLikelyHeart));
fprintf('Complex ICs, N = %d (removed completely).\n',sum(ICsMostLikelyComplex));
fprintf('Total ICs,   N = %d (removed completely).\n',sum(ICsforRemoval));
fprintf('Muscle ICs,  N = %d (keeping >%d Hz).\n',sum(ICsMostLikelyMuscle),muscleFreq);
fprintf('Channel ICs, N = %d (regressed out).\n',sum(ICsMostLikelyChannel));

% 1. Remove selected bad ICs
fprintf('\nRemoving some ICs...\n');
if any(ICsforRemoval)
    artifactComponents(ICsforRemoval,:) = dataICs(ICsforRemoval,:);
end

% 2. Obtain muscle artifact for subtraction by highpass filtering data
fprintf('Filtering muscle ICs...\n');
if any(ICsMostLikelyMuscle)
    [bh, ah] = butter(2, muscleFreq/(EEG.srate/2), 'high');
    artifactComponents(ICsMostLikelyMuscle,:) = do_filteringcore(bh,ah, dataICs(ICsMostLikelyMuscle,:), EEG.event,EEG.srate);
end

% 3. Obtain channel artifact for subtraction by wavelet-thresholding
% -> this did not work well
% if any(ICsMostLikelyChannel)
%     ICsMostLikelyChannel2 = find(ICsMostLikelyChannel);
%
%     for i = 1:length(ICsMostLikelyChannel2)
%         artifactComponents(ICsMostLikelyChannel2(i),:) = wThresholding(dataICs(ICsMostLikelyChannel2(i),:));
%     end
% end

%% Remove artifact and reconstruct data
artifactEEG = EEG.icawinv * artifactComponents;
artifactEEG = reshape(artifactEEG,EEG.nbchan,NTPT,EEG.trials);

% chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% EEGNEW = EEG;
% EEGNEW.data(chaneeg,:) = EEG.data(chaneeg,:) - artifactEEG;
% vis_artifacts(EEGNEW,EEG);

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEG.data(chaneeg,:) = EEG.data(chaneeg,:) - artifactEEG;

%% Regress out channel ICs
fprintf('Regressing out channel ICs...\n');
if any(ICsMostLikelyChannel)
    rawEEG   = EEG.data(chaneeg,:)';
    chanICs  = dataICs(ICsMostLikelyChannel,:)';
    cleanEEG = nt_tsr(rawEEG, chanICs);
    EEG.data(chaneeg,:) = cleanEEG';
end

% EEGNEW = EEG;
% EEGNEW.data(chaneeg,:) = cleanEEG';
% vis_artifacts(EEGNEW,EEG);

%% Log
NBIC = sum(ICsforRemoval | ICsMostLikelyMuscle);

EEG.ALSUTRECHT.ica.ICsforRemoval                       = ICsforRemoval;
EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle2                = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.ICsMostLikelyChannel2               = ICsMostLikelyChannel;

EEG.ALSUTRECHT.ica.numberArtifactICs                   = NBIC;
EEG.ALSUTRECHT.ica.proportionArtifactICs               = NBIC./NICA;
EEG.ALSUTRECHT.ica.numberArtifactICsRemoval            = sum(ICsforRemoval);
EEG.ALSUTRECHT.ica.numberArtifactICsMuscle             = sum(ICsMostLikelyMuscle);
EEG.ALSUTRECHT.ica.numberArtifactICsChannel            = sum(ICsMostLikelyChannel);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'ICA bad component removal\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of bad ICs:    %d\n',EEG.ALSUTRECHT.ica.numberArtifactICs);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of removed IC: %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsRemoval);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of muscle IC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsMuscle);

% Remove (not needed)
EEG.icaact = [];

end