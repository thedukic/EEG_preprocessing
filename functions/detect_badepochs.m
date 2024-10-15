function EEG_cleaned = detect_badepochs(EEG)
% Detect bad epochs using variance and the G-ESD method

% Step 1: Calculate Variance Across Channels for Each Trial
% Get the number of trials
num_trials = size(EEG.data, 3);

% Initialize variance array
variance_per_trial = zeros(1, num_trials);

% Calculate variance across channels for each trial
for trial_idx = 1:num_trials
    trial_data = EEG.data(:, :, trial_idx);
    variance_per_trial(trial_idx) = mean(var(trial_data, 0, 1));
end

% Step 2: Use MATLAB's G-ESD Function to Detect Outliers
max_outliers1 = isoutlier(variance_per_trial, 'gesd');

% Step 3: Compute Temporal Derivative and Repeat
% Compute the temporal derivative of EEG.data along the time axis (2nd dimension)
EEG_derivative = diff(EEG.data, 1, 2);

% Initialize variance array for the derivative data
variance_per_trial_derivative = zeros(1, num_trials);

% Calculate variance across channels for each trial (on the derivative data)
for trial_idx = 1:num_trials
    trial_data_derivative = EEG_derivative(:, :, trial_idx);
    % Note: EEG_derivative will have 1 fewer time point due to diff, but this should be fine
    variance_per_trial_derivative(trial_idx) = mean(var(trial_data_derivative, 0, 1));
end

% Apply G-ESD to the variance of the derivative data
max_outliers2 = isoutlier(variance_per_trial_derivative, 'gesd');

% Step 4: Combine Results and Mark Bad Trials
% Combine bad trials from both the original data and the derivative
combined_bad_trials = max_outliers1 | max_outliers2;

% Display the total number of bad trials
disp(['Total bad trials: ', num2str(sum(combined_bad_trials))]);

% Optional: Remove bad trials from EEG struct if necessary
EEG_cleaned = pop_select(EEG, 'notrial', find(combined_bad_trials));

end