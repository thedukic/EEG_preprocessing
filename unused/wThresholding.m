function data_wt = wThresholding(data)
%
% Function based on: RELAX_wICA_on_ICLabel_artifacts
%
% Input:
%   K = Threshold multiplier for wavelet thresholding.
%       Higher thresh -> Less strict
%   L = Level set for stationary wavelet transform.
%       Higher levels give better freq resolution, but less temp resolution
%   W = Wavelet family to use.
%       Type "wavenames" to see a list of possible wavelets
%
% More info:
% https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2018.00097/full
% Given that the magnitude of artifacts can be far greater than
% that of neurophysiological signals, the component time series
% whose amplitudes are large enough to survive the wavelet-thresholding
% are taken as the artifact timeseries.
%
% Treshold of 0   -> whole IC rejected
% Higher treshold -> Less of the IC is considered as noise
%
% SDukic, January 2025
% =========================================================================
% Params
L = 5;
K = 2;
W = 'coif5'; % coif5 / coif3 / sym8

% Padding
NTPT = length(data);
check_padding_required = mod(NTPT,2^L);
if check_padding_required ~= 0
    padding = zeros(1, (2^L)-check_padding_required);
else
    padding = [];
end

if ~isempty(padding)
    data = [data, padding];
end

% Automatically obtain wavelet enhancement threshold
[wavelet_threshold,threshold_type,~] = ddencmp('den','wv',data);

% Apply stationary wavelet transform to each component to reduce neural contribution to component
wavelet_transform = swt(data,L,W);
wavelet_threshold = K * wavelet_threshold;

% Remove negligible values by applying thresholding
thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold);

% Use inverse wavelet transform to obtained the wavelet transformed component
data_wt = iswt(thresholded_wavelet_transform,W);

% fprintf('\nWavelet tresholding:\n');
% fprintf('Used treshold %1.2f\n',wavelet_threshold);

% Remove padding
if ~isempty(padding)
    data_wt = data_wt(:,1:end-numel(padding));
end

end