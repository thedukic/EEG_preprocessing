function C_zscored = estimate_channelcov(EEG)

% Data
D = EEG.data(1:128,:)';

% Step 1: Compute the correlation matrix C
C = corrcoef(D);

% Step 2: Extract the off-diagonal elements
% Create a mask to ignore the diagonal (which is all 1's)
mask = ~eye(size(C));

% Extract the off-diagonal elements
off_diag_elements = C(mask);

% Step 3: Calculate the mean and standard deviation of the off-diagonal elements
mean_val = mean(off_diag_elements);
std_val  = std(off_diag_elements);

% Step 4: Z-score the matrix
% Z-score all elements of the matrix
C_zscored = (C - mean_val) / std_val;

% Optionally, restore the diagonal to 1's, as we don't z-score the diagonal
C_zscored(logical(eye(size(C)))) = 1;

% Now, C_zscored is the z-scored correlation matrix
figure; imagesc(C_zscored); colorbar; axis square;

end