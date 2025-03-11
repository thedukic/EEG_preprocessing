function EEG = reduce_sparseartifacts(EEG)
% Reduce sparse artifacts
% https://www.sciencedirect.com/science/article/pii/S0165027016000066
%
% -> channel-specific artifats
% -> slow electrode drifts and pops
% -> does not fix EMG, EOG, ECG
% TODO:
% 1. If done on trials, correct the whole trial?
% 2.
%

% Define
Niter  = 2;    % Number of STAR iterations
Nneigh = 12;   % Number of neighouring channels used, keep it 10-15?
Texc   = 3;    % Threshold for excentricity, higher -> looser
Ndeep  = 2;    % Maximum number of channels to fix at each sample
Tpca   = 0.15; % Threshold for discarding weak PCs (percent of the max{PCs} of C of that neigh group of channels)

fprintf('\n================================\n');
fprintf('STAR: Sparse time artifact removal\n');
fprintf('================================\n');

% Check
assert(ismatrix(EEG.data));

% Define neighbours
fprintf('Extracting channel neighbours...\n');
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
roiNeighbours = find_neighbours(EEG.chanlocs(chaneeg),100);
roiNeighbours = roiNeighbours(:,1:(Nneigh+1));

% The first column should be the electrode itself
% -> Remove that column
NCHN = sum(chaneeg);
maskNeigh = NaN(NCHN,1);
for i = 1:NCHN
    maskNeigh(i) = find(roiNeighbours(i,:) == i);
end

assert(~any(isnan(maskNeigh)));
assert(all(maskNeigh == 1));
roiNeighbours(:,1) = [];

% STAR
fprintf('Applying %d iterations.\n',Niter);
x = double(EEG.data(chaneeg,:))';

for i = 1:Niter
    fprintf('Iteration: %d\n',i);
    [x, w, ww] = nt_star(x,Texc,roiNeighbours,Ndeep,Tpca);
end

% % Check
% EEGNEW = EEG;
% EEGNEW.data(chaneeg,:) = x';
% vis_artifacts(EEGNEW,EEG);

% Store
EEG.data(chaneeg,:) = x';

end