function y = robust_zscore(dataToStandardize, referenceData)
%ROBUST_ZSCORE Performs robust Z-score standardization.
%
%   Y = ROBUST_ZSCORE(X) standardizes the data in X using robust measures
%   of central tendency (median) and spread (Median Absolute Deviation, MAD).
%   The MAD is scaled by 1.4826 to be a consistent estimator of the
%   standard deviation for normally distributed data.
%   If X is a matrix, the operation is performed column-wise.
%
%   Y = ROBUST_ZSCORE(X, REF) standardizes the data in X using the median
%   and MAD calculated from the reference data REF. This is useful for
%   standardizing new data points (X) relative to a control or baseline
%   dataset (REF) without letting outliers in X influence the scaling.
%   If REF is a matrix, its median and MAD are calculated column-wise.
%
%   References:
%   [] W.A. Stahel. Robuste Schatzungen: infinitesimale Optimalit¨at und
%      Schatzungen von Kovarianzmatrizen. PhD thesis, ETH Zurich, 1981.
%   [] D.L. Donoho. Breakdown properties of multivariate location estimators.
%      Qualifying paper, Harvard University, Boston, 1982.
%
%   Inputs:
%     dataToStandardize - Numeric array (vector or matrix) to be standardized.
%     referenceData     - (Optional) Numeric array (vector or matrix)
%                         from which median and MAD are calculated.
%                         If omitted, dataToStandardize is used as its own reference.
%
%   Output:
%     y - Standardized data, same size as dataToStandardize.
%         Columns with zero MAD will result in NaN for that column's output.

% Input validation
if ~isnumeric(dataToStandardize) || ~ismatrix(dataToStandardize)
    error('robust_zscore:InvalidInput', 'Input ''dataToStandardize'' must be a numeric vector or matrix.');
end
if isempty(dataToStandardize)
    y = []; % Return empty if input is empty
    return;
end

% Determine if a reference dataset is provided
useReference = (nargin == 2);

if useReference
    if ~isnumeric(referenceData) || ~ismatrix(referenceData)
        error('robust_zscore:InvalidReferenceInput', 'Input ''referenceData'' must be a numeric vector or matrix.');
    end
    if isempty(referenceData)
        error('robust_zscore:EmptyReference', 'Input ''referenceData'' cannot be empty when provided.');
    end

    % Calculate robust statistics from the reference data (column-wise)
    robust_median = median(referenceData, 'omitnan');
    raw_mad = mad(referenceData, 1, 'omitnan'); % 1 for median-based MAD
else
    % Calculate robust statistics from the data to be standardized (column-wise)
    robust_median = median(dataToStandardize, 'omitnan');
    raw_mad = mad(dataToStandardize, 1, 'omitnan'); % 1 for median-based MAD
end

% Scaling factor for MAD to estimate standard deviation for normal data
scaling_factor_for_std = 1.4826;
robust_std_estimate = scaling_factor_for_std * raw_mad;

% Handle cases where robust_std_estimate is zero (i.e., all values in a column are identical)
% In such cases, standardization is not meaningful, so return NaN for that column.
% This is generally preferred over returning Inf.
zero_mad_cols = (robust_std_estimate == 0);
if any(zero_mad_cols)
    warning('robust_zscore:ZeroMAD', ...
        'One or more columns have a zero Median Absolute Deviation. Corresponding output will be NaN.');
end

% Perform the robust z-scoring
y = (dataToStandardize - robust_median) ./ robust_std_estimate;

% For columns where MAD was zero, set output to NaN
y(:, zero_mad_cols) = NaN;

end


% function y = robust_zscore(x1,x2)
% % [] W.A. Stahel. Robuste Schatzungen: infinitesimale Optimalit¨at und Schatzungen von Kovarianzmatrizen. PhD thesis, ETH Zurich, 1981.
% % [] D.L. Donoho. Breakdown properties of multivariate location estimators. Qualifying paper, Harvard University, Boston, 1982.
% 
% % x1, x2 are either vectirs or matrices
% % If a matrix (eg Nsubj x Neeg), then oparates along the columns
% 
% if nargin == 1
%     % disp('Robust Z-score using patient data.');
%     y = (x1-median(x1)) ./ (1.4826*mad(x1,1));     % median absolute deviation
%     % y = (x1-median(x1)) ./(1.253314*mad(x1,0));  % mean absolute deviation
% else
%     % disp('Robust Z-score using control data.');
%     y = (x1-median(x2)) ./ (1.4826*mad(x2,1));      % median absolute deviation
%     % y = (x1-median(x2)) ./ (1.253314*mad(x2,0));  % mean absolute deviation
% end
% 
% end