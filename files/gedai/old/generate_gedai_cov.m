function generate_gedai_cov
% Script to Calculate Reference Covariance Matrices from Leadfields
%
% This script iterates through all .mat files in a specified directory,
% loads the 'leadfield' structure from each file, extracts the leadfield
% matrix (L) for inside vertices, centres the matrix, normalises it to
% ensure equal subject contribution, and calculates the covariance matrix
% (L * L') for reference purposes.

% -------------------------------------------------------------------------
% 1. Configuration
% -------------------------------------------------------------------------
folderPath_individual = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\individual\leadfields';
folderPath_template   = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\template\eeg\BEM';
folderPath_save       = 'C:\DATA\MATLAB\myCodes\preprocessing\files\gedai';

% -------------------------------------------------------------------------
% 2. File Finding and Setup
% -------------------------------------------------------------------------
% Find all .mat files in the specified folder
fileList = dir(fullfile(folderPath_individual, '*.mat'));
if isempty(fileList)
    error('No .mat files found in the specified directory: %s', folderPath_individual);
end
numFiles = length(fileList);
fprintf('Found %d leadfield files to process.\n', numFiles);

% Initialise the 3D array for storing the covariance matrices.
isFirstPass = true;
num_channel = 128;

% -------------------------------------------------------------------------
% 3. Processing Loop
% -------------------------------------------------------------------------
for i = 1:numFiles
    fileName = fileList(i).name;
    fullFilePath = fullfile(folderPath_individual, fileName);
    fprintf('Processing file %d/%d: %s\n', i, numFiles, fileName);
    try
        % Load the structure named 'leadfield' from the .mat file
        leadfield = load(fullFilePath, 'leadfield');
        leadfield = leadfield.leadfield;

        % (a) Identify lead sources that are 'inside' the brain surface
        indx = find(leadfield.inside);

        % (b) Extract and concatenate the leadfield vectors for these sources.
        L = cell2mat(leadfield.leadfield(indx).');

        % (c) Centre the leadfield matrix L by removing the mean of each column (dipole orientation)
        % L = L - mean(L, 1);
        L_mean = sum(L, 1) / (num_channel + 1);
        L = L - L_mean;

        % (d) Normalise the leadfield matrix by its Frobenius norm
        % EXPLANATION: This step ensures that ||L||_F = 1. By doing this before
        % calculating the covariance, we guarantee that every subject contributes
        % exactly the same amount of 'power' to the final average reference matrix.
        % This prevents subjects with a higher number of internal vertices or
        % different electrode geometries from dominating the reference spatial pattern.
        L = L / norm(L, 'fro');

        % (e) Calculate the covariance matrix: L * L'
        % This results in a N_sensors x N_sensors matrix (Covariance across sensors).
        covMat = L * L';

        % -----------------------------------------------------------------
        % 4. Storage and Pre-allocation
        % -----------------------------------------------------------------
        if isFirstPass
            % Pre-allocate the covRef array based on the size of the first matrix
            covRef_all = NaN(num_channel, num_channel, numFiles+1);
            isFirstPass = false;
            assert(size(L, 1) == num_channel);
        end
        % Store the resulting covariance matrix in the 3D array
        covRef_all(:, :, i) = covMat;
    catch ME
        fprintf('An error occurred while processing %s: %s\n', fileName, ME.message);
        % Display the error identifier for debugging
        disp(ME.identifier);
        continue;
    end
end

% -------------------------------------------------------------------------
% 5. Add the template too
% -------------------------------------------------------------------------
fullFilePath = fullfile(folderPath_template, 'leadfield_vrtx.mat');
load(fullFilePath, 'leadfield');

indx = find(leadfield.inside);
L = cell2mat(leadfield.leadfield(indx).');
% L = L - mean(L, 1);
L_mean = sum(L, 1) / (num_channel + 1);
L = L - L_mean;

% Normalise the template leadfield as well
% EXPLANATION: Ensures the template matrix is on the exact same scale
% as the individual subjects before it is added to the array.
L = L / norm(L, 'fro');

covMat = L * L';
covRef_all(:, :, numFiles+1) = covMat;

% -------------------------------------------------------------------------
% 6. Final Output
% -------------------------------------------------------------------------
fullFilePath = fullfile(folderPath_save, 'leadfield_template.mat');
save(fullFilePath, 'L');

% Remove any trailing zero-slices if some files were skipped
fprintf('\nProcessing complete. Successfully calculated %d covariance matrices.\n', size(covRef_all, 3));
fprintf('The final covariance reference matrix array is stored in the variable ''covRef'' (Size: %dx%dx%d).\n', size(covRef_all));

% Save
covRef = mean(covRef_all, 3);
% covRef = average_car_covariances(covRef_all);
fullFilePath = fullfile(folderPath_save, 'covRef_template+individual.mat');
save(fullFilePath, "covRef");

% Save
covRef = covRef_all(:, :, end);
fullFilePath = fullfile(folderPath_save, 'covRef_template.mat');
save(fullFilePath, "covRef");

% -------------------------------------------------------------------------
% 6. Final Output
% -------------------------------------------------------------------------
% Assuming you have calculated both averages
covMeanEuclidean = mean(covRef_all, 3);
covMeanLogEuc = average_car_covariances(covRef_all);

% Extract eigenvalues and sort them in descending order
eig_Euc = sort(eig(covMeanEuclidean), 'descend');
eig_Log = sort(eig(covMeanLogEuc), 'descend');

% Plot the spectrum (ignoring the 128th zero eigenvalue)
figure; hold on;
% Using a logarithmic Y-axis because the first few eigenvalues are massive
semilogy(1:127, eig_Euc(1:127), 'LineWidth', 2, 'DisplayName', 'Euclidean (Swollen)');
semilogy(1:127, eig_Log(1:127), 'LineWidth', 2, 'DisplayName', 'Log-Euclidean (True)');

title('Eigenvalue Spectrum Comparison');
xlabel('Component Number');
ylabel('Variance (Log Scale)');
legend; grid on;

% Extract the dominant spatial pattern (First Eigenvector)
[V_Euc, D_Euc] = eigs(covMeanEuclidean, 1);
[V_Log, D_Log] = eigs(covMeanLogEuc, 1);

% Ensure both vectors are flipped in the same direction for visual comparison
if sign(max(V_Euc)) ~= sign(max(V_Log))
    V_Log = -V_Log;
end

% Plot using your existing topoplot function
figure; tiledlayout(1,2, 'TileSpacing', 'compact');
mytopoplot(V_Euc, [], 'Euclidean Average', nexttile);
mytopoplot(V_Log, [], 'Log-Euclidean Average', nexttile);

end

function covMeanLE = average_car_covariances(covMatrix3D)
% Averages Common Average Referenced (CAR) covariance matrices
% using Subspace Projection to avoid non-positive eigenvalues.

[nChan, ~, nSubj] = size(covMatrix3D);

% Step 1: Temporarily drop the last channel to make it full rank (nChan - 1)
% Does not matter which one is removed!
subCovMatrix = covMatrix3D(1:end-1, 1:end-1, :);
sumLogMat = zeros(nChan-1, nChan-1);

% Step 2: Log-Euclidean average in the reduced SPD subspace
for i = 1:nSubj
    currentCov = subCovMatrix(:,:,i);
    currentCov = (currentCov + currentCov') / 2; % Ensure symmetry

    logMat = logm(currentCov);
    logMat = (logMat + logMat') / 2;

    sumLogMat = sumLogMat + logMat;
end

meanLogMat = sumLogMat / nSubj;
subCovMean = expm(meanLogMat);
subCovMean = (subCovMean + subCovMean') / 2;

% Step 3 & 4: Reconstruct using the Interpolation Matrix (T)
% T is size [nChan x nChan-1]
T = [eye(nChan - 1); -ones(1, nChan - 1)];

% Project the full-rank subspace back to the CAR constraint
covMeanLE = T * subCovMean * T';

% Enforce final symmetry
covMeanLE = (covMeanLE + covMeanLE') / 2;

end