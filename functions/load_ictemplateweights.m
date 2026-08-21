function templates_ica = load_ictemplateweights(DATA)

fprintf('--------------------------------\n');
fprintf('Loading ICA templates\n');
fprintf('--------------------------------\n');

% Load individual tempaltes
fileName = [DATA.ALSUTRECHT.subject.id '_' DATA.ALSUTRECHT.subject.visit '_' DATA.ALSUTRECHT.subject.task '_wartifacts.mat'];
load(fullfile(DATA.ALSUTRECHT.subject.data, fileName), "wArtifacts");

% Load group tempaltes
load('wBlink.mat', 'Blinkweights');
load('wSaccade.mat', 'Saccadeweights');
load('wHeart.mat', 'Heartweights');

% Organise
templates_ica.Blinkweights0    = Blinkweights(1:128);
templates_ica.Blinkweights1    = wArtifacts.Blinkweights1(1:128);

templates_ica.Saccadeweights0  = Saccadeweights(1:128);
templates_ica.Saccadeweights1L = wArtifacts.Saccadeweights1L(1:128);
templates_ica.Saccadeweights1R = wArtifacts.Saccadeweights1R(1:128);

templates_ica.Heartweights0    = Heartweights(1:128, :);
templates_ica.Heartweights1    = wArtifacts.Heartweights1(1:128);
templates_ica.Heartweights2    = wArtifacts.Heartweights2(1:128);

fprintf('Done!\n');

% % --- Blink Weights ---
% if wArtifacts.Blinkcorr < 0.6
%     load('wBlink.mat','Blinkweights');
%     fprintf('Blink weights: Low correlation (%.2f). Loading template Blinkweights.\n', wArtifacts.Blinkcorr);
% else
%     Blinkweights = wArtifacts.Blinkweights(1:128);
%     fprintf('Blink weights: High correlation (%.2f). Using individual Blinkweights.\n', wArtifacts.Blinkcorr);
% end
%
% % --- Saccade Weights ---
% if wArtifacts.Saccadecorr1 < 0.6
%     load('wSaccade.mat','Saccadeweights');
%     fprintf('Saccade weights 1: Low correlation (%.2f). Loading template Saccadeweights.\n', wArtifacts.Saccadecorr1);
%     Saccadeweights1 = Saccadeweights;
% else
%     Saccadeweights1 = wArtifacts.Saccadeweights1(1:128);
%     fprintf('Saccade weights 1: High correlation (%.2f). Using individual Saccadeweights.\n', wArtifacts.Saccadecorr1);
% end
% if wArtifacts.Saccadecorr2 < 0.6
%     load('wSaccade.mat','Saccadeweights');
%     fprintf('Saccade weights 2: Low correlation (%.2f). Loading template Saccadeweights.\n', wArtifacts.Saccadecorr2);
%     Saccadeweights2 = Saccadeweights;
% else
%     Saccadeweights2 = wArtifacts.Saccadeweights2(1:128);
%     fprintf('Saccade weights 2: High correlation (%.2f). Using individual Saccadeweights.\n', wArtifacts.Saccadecorr2);
% end
%
% % --- Heart Weights ---
% if wArtifacts.Heartcorr < 0.6
%     load('wHeart.mat','Heartweights');
%     fprintf('Heart weights: Low correlation (%.2f). Loading template Heartweights.\n', wArtifacts.Heartcorr);
% else
%     Heartweights = wArtifacts.Heartweights(1:128);
%     fprintf('Heart weights: High correlation (%.2f). Using individual Heartweights.\n', wArtifacts.Heartcorr);
% end

end