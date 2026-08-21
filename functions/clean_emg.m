function EMG = clean_emg(EMG)

% -------------------------------------------------------------------------
% Median filter
% -------------------------------------------------------------------------
% % Ensure EMG structure has data and srate
% if ~isfield(EMG, 'data') || ~isfield(EMG, 'srate')
%     error('EMG structure must contain fields ''data'' and ''srate''.');
% end
% 
% fprintf('\n================================\n');
% fprintf('Removing spikes from EMG data\n');
% fprintf('================================\n');
% 
% % --- 1. Define Time-Based Window ---
% % Target duration for spike removal (e.g., 5 milliseconds or 0.005 seconds)
% T_target = 0.005;
% 
% % Calculate required window size (N) based on sampling frequency
% % Round to the nearest odd integer to maintain temporal consistency across srate
% window_size_calc = round(T_target * EMG.srate);
% 
% % Ensure the window size is an odd number (median filter convention)
% if mod(window_size_calc, 2) == 0
%     window_size = window_size_calc + 1;
% else
%     window_size = window_size_calc;
% end
% 
% % Ensure minimum window size is 3
% window_size = max(window_size, 3);
% 
% % --- 2. Filtering ---
% % Determine the time dimension (assuming time is the dimension with more samples)
% if size(EMG.data, 1) > size(EMG.data, 2)
%     time_dim = 1; % Time is in rows
% else
%     time_dim = 2; % Time is in columns (most common for multi-channel data)
% end
% 
% fprintf('-> Sampling Frequency (srate): %.0f Hz\n', EMG.srate);
% fprintf('-> Calculated Window Size (N): %d samples (%.2f ms)\n', window_size, window_size * 1000 / EMG.srate);
% 
% % Apply the median filter along the time dimension
% % The time_dim argument makes the function robust to data orientation.
% EMG.data = medfilt1(EMG.data, window_size, [], time_dim);
% 
% fprintf('Done! Spikes removed.\n');

% -------------------------------------------------------------------------
% STAR
% -------------------------------------------------------------------------
EMG = do_star(EMG, 'emg');

end

