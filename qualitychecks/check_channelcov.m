function EEG = check_channelcov(EEG)

doZscore = false;

% fprintf('Estimating a correlation matrix from the cleaned data...\n');

% Select channels
eegchan = strcmp({EEG.chanlocs.type},'EEG');
% extchan = ismember({EEG.chanlocs.labels},{'VEOG','HEOG','LM','RM'});
extchan = ismember({EEG.chanlocs.labels},{'VEOG','HEOG','ECG'});
selchan = eegchan | extchan;

% Data
D = EEG.data(selchan,:,:);
D = permute(D,[2 1 3]);

% Compute the correlation matrix C
NTRL = EEG.trials;
NCHN = sum(selchan);
% NLEN = EEG.pnts-1;

C = NaN(NCHN,NCHN,NTRL);
for i = 1:NTRL
    % datTmp = D(:,:,i);
    % datTmp = datTmp-mean(datTmp,2);
    % C(:,:,i) = (datTmp*datTmp')./NLEN;
    C(:,:,i) = corr(D(:,:,i));
end

% Average
C = mean(C,3);

% Zscoring is not needed here and should be avoided?
if doZscore
    % Extract the off-diagonal elements
    % Create a mask to ignore the diagonal (which is all 1's)
    mask = ~eye(size(C));

    % Extract the off-diagonal elements
    off_diag_elements = C(mask);

    % Calculate the mean and standard deviation of the off-diagonal elements
    mean_val = mean(off_diag_elements);
    std_val  = std(off_diag_elements);

    % Z-score the matrix
    C = (C - mean_val) / std_val;
end

% Ensure the diagonal has 1's
C(logical(eye(size(C)))) = 1;

% Log
EEG.ALSUTRECHT.chanCorr = C;

% % Plot
% figure; imagesc(C); colorbar; axis square;

end