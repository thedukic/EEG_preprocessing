function EEG3 = do_relica(EEG,subject)
% https://www.sciencedirect.com/science/article/pii/S1053811914007526

rankICA = 20;
bootICA = 50;

% -------------------------------------------------------------------------
% 1. Reshape the data to 2D (channels by total timepoints)
% 2. Mean-centre the data (required for accurate Principal Component Analysis)
assert(ismatrix(EEG.data));
EEG.data = double(EEG.data);
EEG.data = EEG.data - mean(EEG.data, 2);

% 3. Calculate the Principal Component Analysis
% 'coeff' contains the spatial weights, 'score' contains the time courses
[coeff, score, ~] = pca(EEG.data');

% 4. Select only the first 'rankICA' components
pca_weights = coeff(:, 1:rankICA);

% 5. Project the original data into the reduced component space
% The new data size will be (rankICA by total timepoints)
reduced_data = pca_weights' * EEG.data;

% 6. Create the new EEG2 structure
EEG2 = EEG;

% 7. Replace the data and update the channel dimensions
EEG2.data = reshape(reduced_data, rankICA, EEG.pnts, EEG.trials);
EEG2.nbchan = rankICA;

% 8. Remove physical channel locations as the "channels" are now abstract components
EEG2.chanlocs = [];

% -------------------------------------------------------------------------
% 9. Run RELICA on the reduced dataset
EEG2 = relica(...
    EEG2, bootICA, 'beamica', 'point', subject.preproc, 'local', ...
    'icaopt', {'extended', 0});

% -------------------------------------------------------------------------
EEG2 = relica_plots(EEG2, 'cluster');

% -------------------------------------------------------------------------
% Select the centroid
EEG2.icawinv     = EEG2.etc.RELICA.A_centroid;
EEG2.icaweights  = EEG2.etc.RELICA.W_centroid;
EEG2.icasphere   = eye(EEG2.nbchan);
EEG2.icachansind = 1:EEG2.nbchan;
EEG2 = eeg_checkset(EEG2);

% -------------------------------------------------------------------------
% Calculate the final unmixing matrix for the physical scalp electrodes
final_weights = EEG2.icaweights * EEG2.icasphere * pca_weights';

% Create a new Electroencephalography structure to preserve the original data
EEG3 = EEG;
EEG3.etc.RELICA = EEG2.etc.RELICA;

% Apply the back-projected weights to the new structure
EEG3.icaweights = final_weights;
EEG3.icasphere  = eye(size(final_weights, 2)); % Set sphere to an identity matrix
EEG3.icawinv    = pinv(final_weights);         % Calculate the topographical spatial projections

% Define which physical channels were used in the analysis
% Assuming you used all channels. If you excluded some (like Electrooculography channels), update this vector.
EEG3.icachansind = 1:EEG3.nbchan;

% Run the internal validation function to ensure structure consistency
EEG3 = eeg_checkset(EEG3);

% -------------------------------------------------------------------------
% pop_topoplot(EEG3,0);

end