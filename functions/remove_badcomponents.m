function EEG = remove_badcomponents(EEG)
%
%
% SDukic, November 2024
% =========================================================================

% Muscle ICs will be filtered such that they have only freq above this one
muscleFreq = 25; % [Hz]

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

% Testing
fprintf('Muscle ICs above the 25th IC are removed!\n');
ICsMostLikelyMuscle(26:end) = false;

% *Remove completely ICs
ICsforRemoval = false(NICA,1);
ICsforRemoval(ICsMostLikelyEye | ICsMostLikelyHeart  | ICsMostLikelyComplex) = true;

% Report
fprintf('Eye ICs,     N = %d (removed completely).\n',sum(ICsMostLikelyEye));
fprintf('Heart ICs,   N = %d (removed completely).\n',sum(ICsMostLikelyHeart));
fprintf('Complex ICs, N = %d (removed completely).\n',sum(ICsMostLikelyComplex));
fprintf('Total ICs,   N = %d (removed completely).\n',sum(ICsforRemoval));
fprintf('Muscle ICs,  N = %d (keeping >%d Hz).\n',sum(ICsMostLikelyMuscle),muscleFreq);

% 1. Remove selected bad ICs
if any(ICsforRemoval)
    artifactComponents(ICsforRemoval,:) = dataICs(ICsforRemoval,:);
end

% 2. Obtain muscle artifact for subtraction by highpass filtering data
if any(ICsMostLikelyMuscle)
    [bh, ah] = butter(2, muscleFreq/(EEG.srate/2), 'high');
    artifactComponents(ICsMostLikelyMuscle,:) = do_filteringcore(bh,ah,dataICs(ICsMostLikelyMuscle,:),EEG.event,EEG.srate);
end

%% Remove artifact and reconstruct data:
artifactEEG = EEG.icawinv * artifactComponents;
artifactEEG = reshape(artifactEEG,EEG.nbchan,NTPT,EEG.trials);

% chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% EEGNEW = EEG;
% EEGNEW.data(chaneeg,:,:) = EEG.data(chaneeg,:,:) - artifactEEG;
% vis_artifacts(EEGNEW,EEG);

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:) - artifactEEG;

%% Log
NBIC = sum(ICsforRemoval | ICsMostLikelyMuscle);

EEG.ALSUTRECHT.ica.ICsforRemoval                       = ICsforRemoval;
EEG.ALSUTRECHT.ica.numberArtifactICs                   = NBIC;
EEG.ALSUTRECHT.ica.proportionArtifactICs               = NBIC./NICA;
EEG.ALSUTRECHT.ica.numberArtifactICsRemoval            = sum(ICsforRemoval);
EEG.ALSUTRECHT.ica.numberArtifactICsMuscle             = sum(ICsMostLikelyMuscle);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA cleaning\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of bad ICs:    %d\n',EEG.ALSUTRECHT.ica.numberArtifactICs);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of removed IC: %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsRemoval);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of muscle IC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsMuscle);

% Not needed
EEG.icaact = [];

end