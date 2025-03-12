function NICA = select_nICs(EEG,NICA)
%
% Estimating the number of PCs that will be used for ICA
% See: https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline#What_is_the_minimum_length_of_data_to_perform_ICA.3F_.2807.2F04.2F2022_adde
%
% EEG available:
% MMN  3*7     ~ 21 min
% SART 3*5     ~ 15 min
% RS   2x3x2   ~ 12 min
% MT   7+3+7   ~ 17 min
%
% Minimum EEG needed:
% 128 elecs ^2*30 /256/60 ~ 32  min
% 70 PCs    ^2*30 /256/60 ~ 9.5 min
% 50 PCs    ^2*30 /256/60 ~ 5.0 min
%
% -> estimate 70 IC unless to little data, then estimate 50 ICs
% -> the latter was added to account for DUB EO data (~ 6 min)

fprintf('\n--------------------------------\n');
fprintf('Selecting the number of PCs for ICA\n');
fprintf('--------------------------------\n');

dataLength = prod(size(EEG.data,[2 3]));
NICAmax = estim_optimalN(dataLength);

if dataLength > estim_dataneeded(NICA(1))
    NICA = NICA(1);
elseif dataLength > estim_dataneeded(NICA(2))
    NICA = NICA(2);
else
    warning('The data seems to be very short.');
    NICA = NICA(2) - 10;
end

% Report
fprintf('Number of PCs for ICA: %d\n', NICA);
fprintf('Max number of PCs for ICA: %d\n', NICAmax);

end

% =========================================================================
% Helper function
% =========================================================================
function L = estim_dataneeded(NICA)
L = 30 * (NICA^2);
end

function NICA = estim_optimalN(L)
NICA = round(sqrt((L/30)));
end