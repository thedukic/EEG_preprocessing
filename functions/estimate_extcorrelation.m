function R = estimate_extcorrelation(DATA, maskArtifact, typeArtifact)
% ESTIMATE_EXTCORRELATION computes the Spearman correlation matrix between
% all EEG channels and the specified external artifact channels contained in an EEGLAB structure.
%
% This function is typically used to derive weights (correlation values)
% that inform subsequent artifact correction methods (e.g., AAR, regression)
% by identifying which EEG channels are most contaminated by external
% physiological artifacts (EOG, ECG, etc.).
%
% Inputs:
%   DATA (struct): An EEGLAB-format structure containing the data and channel information.
%                  Expected fields:
%                  - DATA.data (Full data matrix: N_Channels x N_Time_Points)
%                  - DATA.chanlocs (Channel locations structure, used to identify EEG/Artifact channels)
%   typeArtifact (string or cell array of strings): The label(s) of the external artifact
%                 channel(s) to correlate with the EEG data (e.g., 'VEOG' or {'VEOG', 'HEOG', 'ECG'}).
%
% Output:
%   R (matrix): Spearman correlation values. Dimensions are
%               (N_EEG_Channels x N_Artifact_Channels), where:
%               R(i, j) is the correlation between EEG channel 'i' and Artifact channel 'j'.
%

% --- 1. Data Extraction and Preparation from EEGLAB structure ---
% Identify EEG channels (type 'EEG').
eeg_idx = strcmp({DATA.chanlocs.type}, 'EEG');

% Identify external artifact channels by matching their labels to typeArtifact.
% We use ismember to correctly handle when typeArtifact is a single string or a cell array of strings.
ext_idx = ismember({DATA.chanlocs.labels}, upper(typeArtifact));

% Input Validation
if sum(eeg_idx) == 0
    error('Could not find any channels marked as type ''EEG'' in DATA.chanlocs.');
end
if sum(ext_idx) == 0
    % If there are no external channels, the correlation is moot.
    warning('No external artifact channels found (i.e., no channels not marked as type ''EEG''). Returning empty matrix.');
    R = [];
    return;
end

% Take all data if the mask is empty
if isempty(maskArtifact)
    maskArtifact = true(1, DATA.pnts);
end

% Separate the data matrices
brainData = DATA.data(eeg_idx, maskArtifact);
extData = DATA.data(ext_idx, maskArtifact);

% Check if the time point dimensions match
if size(brainData, 2) ~= size(extData, 2)
    error('The number of time points (columns) in brainData and extData must be equal.');
end

% --- 2. Correlation Calculation ---
fprintf('Calculating Spearman correlation matrix for %s...\n', typeArtifact);

% The built-in MATLAB 'corr' function calculates correlation between columns.
% To correlate channels (which are rows in our input matrices) across time,
% we must transpose the matrices before calculating.
% R will have size (N_EEG_Channels x N_Artifact_Channels)

R = corr(brainData', extData', 'Type', 'Spearman');

% --- 3. Display Results (Optional) ---
% % --- Display Topoplot (Requires EEGLAB Toolbox) ---
% % The output R is a rectangular matrix (EEG x Artifacts), so a standard
% % topoplot is not suitable unless you plot the correlation vector for a
% % single artifact channel (e.g., R(:, 1) for VEOG).
% % The line below from the original code would only work if the EEGLAB
% % toolbox and a compatible environment is loaded, and typically 'R' should be
% % a vector for topoplot. We leave it commented out as per your previous file.
% % figure; mytopoplot(R,[],'R',nexttile);

end