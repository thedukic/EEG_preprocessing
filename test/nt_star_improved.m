function [x, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH)
% [y, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH)
%
% Sensor noise suppression via Sparse Time Artifact Removal (STAR).
%
% y:          Denoised data
% w:          0 for parts needing fixing, 1 elsewhere (time*1)
% ww:         0 for parts needing fixing, 1 elsewhere (time*chans)
%
% x:          Data to denoise (time*chans or time*chans*trials)
% thresh:     Threshold for excentricity measure (default: 1)
% closest:    Indices of closest channels per row (default: all)
% depth:      Maximum number of channels to fix per sample (default: 1)
% PCA_THRESH: Threshold for discarding weak PCs (default: 1e-15)
% NSMOOTH:    Smoothing window in samples for excentricity (default: 10)

MINPROP = 0.3;  % Minimum proportion of artifact-free data at first iteration
NITER   = 3;    % Iterations to refine c0
VERBOSE = 1;

if nargin < 1, error('More inputs needed.'); end
if nargin < 2 || isempty(thresh), thresh = 1; end
if nargin < 3, closest = []; end
if nargin < 4 || isempty(depth), depth = 1; end
if nargin < 5 || isempty(PCA_THRESH), PCA_THRESH = 1e-15; end
if nargin < 6 || isempty(NSMOOTH), NSMOOTH = 10; end

if ~isempty(closest) && size(closest, 1) ~= size(x, 2)
    error('closest array should have as many rows as channels of x');
end
if any(closest < 0, 'all'), error('Invalid channel indices.'); end

if nargout == 0
    [y, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH);
    disp([mean(w(:)), mean(ww(:))]);
    figure(1); clf;
    subplot(3, 1, 1); plot(x); title('Original');
    subplot(3, 1, 2); plot(y); title('Cleaned');
    subplot(3, 1, 3); plot(w, '.-'); ylim([0 1.1]); title('Clean Mask w');
    clear x w ww;
    return;
end

[nsample, nchan, ~] = size(x);
x = nt_unfold(x);

p0 = nt_wpwr(x);
mn = mean(x, 1);
x  = nt_demean(x);
nn = sqrt(mean(x.^2, 1));
nn(nn == 0) = 1; % Guard against division by zero on dead channels
x  = x ./ nn;
p00 = nt_wpwr(x);

% Flag and isolate flat / NaN channels
iNan  = find(all(isnan(x), 1));
iZero = find(all(x == 0, 1));
x(:, iNan)  = randn(size(x, 1), numel(iNan));
x(:, iZero) = randn(size(x, 1), numel(iZero));

x  = nt_demean(x);
c0 = nt_cov(x);

% Parse neighbour list for each channel
oChanList = cell(nchan, 1);
for iChan = 1:nchan
    if ~isempty(closest)
        oC = closest(iChan, :);
    else
        oC = setdiff(1:nchan, iChan);
    end
    oC(oC > nchan | oC <= 0) = [];
    oChanList{iChan} = oC;
end

% -------------------------------------------------------------------------
% Phase 1: Iterative clean mask and covariance estimation
% -------------------------------------------------------------------------
iIter = NITER;
while iIter > 0
    w = ones(size(x, 1), 1);
    
    for iChan = 1:nchan
        oC = oChanList{iChan};
        if isempty(oC), continue; end
        
        c_sub = c0(oC, oC);
        [topcs, eigenvalues] = nt_pcarot(c_sub);
        keep = (eigenvalues / max(eigenvalues)) > PCA_THRESH;
        topcs = topcs(:, keep);
        
        % Projection vector: x(:, oC) * P_chan predicts x(:, iChan)
        b = (c0(iChan, oC) * topcs) / (topcs' * c_sub * topcs);
        P_chan = topcs * b'; 
        
        y  = x(:, oC) * P_chan;
        dx = abs(y - x(:, iChan)) + eps;
        
        % Excentricity relative to current clean samples
        d = dx / mean(dx(w > 0));
        if NSMOOTH > 0
            d = movmean(d, NSMOOTH, 1);
        end
        
        d = d / thresh;
        w = min(w, double(d < 1));
    end
    
    prop = mean(w);
    if VERBOSE > 0, disp(['Proportion artifact free: ', num2str(prop)]); end
    
    if iIter == NITER && prop < MINPROP
        thresh = 1.1 * thresh;
        if VERBOSE > 0, disp(['Warning: nt_star increasing threshold to ', num2str(thresh)]); end
        w = ones(size(w));
    else
        iIter = iIter - 1;
    end
    
    x  = nt_demean(x, w);
    c0 = nt_cov(x, [], w);
end

% -------------------------------------------------------------------------
% Phase 2: Excentricity ranking and local channel repair
% -------------------------------------------------------------------------
std_clean = sqrt(mean(x(w > 0, :).^2, 1));
std_clean(std_clean == 0) = 1;
xx = abs(x) ./ std_clean;

if NSMOOTH > 0
    xx = movmean(xx, NSMOOTH, 1);
end

% Precompute projection vectors for all channels (computed once from stable c0)
projWeights = cell(nchan, 1);
for iChan = 1:nchan
    oC = oChanList{iChan};
    if isempty(oC), continue; end
    
    c_sub = c0(oC, oC);
    [topcs, eigenvalues] = nt_pcarot(c_sub);
    keep = (eigenvalues / max(eigenvalues)) > PCA_THRESH;
    topcs = topcs(:, keep);
    
    b = (c0(iChan, oC) * topcs) / (topcs' * c_sub * topcs);
    projWeights{iChan} = topcs * b';
end

% Identify worst channels at contaminated points
depth = min(depth, nchan - 1);
[~, rank] = sort(xx, 2, 'descend');
rank(w > 0, :) = 0; % Do not touch samples confirmed artifact-free

ww = ones(size(x));
for iDepth = 1:depth
    iFixed = nchan;
    nFixed = 0;
    
    for iChan = 1:nchan
        bad_samples = find(rank(:, iDepth) == iChan);
        
        if iDepth > 1 && ~isempty(bad_samples)
            bad_samples(xx(bad_samples, iChan) < thresh) = [];
        end
        
        if isempty(bad_samples)
            iFixed = iFixed - 1;
            continue;
        end
        
        nFixed = nFixed + numel(bad_samples);
        ww(bad_samples, iChan) = 0;
        
        oC = oChanList{iChan};
        if ~isempty(oC) && ~isempty(projWeights{iChan})
            % Fix bad samples using precomputed spatial projection
            x(bad_samples, iChan) = x(bad_samples, oC) * projWeights{iChan};
        end
    end
    
    if VERBOSE > 0
        disp(['Depth: ', num2str(iDepth), ...
              ', channels active: ', num2str(iFixed), ...
              ', samples repaired: ', num2str(nFixed), ...
              ', power ratio: ', num2str(nt_wpwr(x) / p00)]);
    end
end

% Restore channel scaling, offsets, and original structure
x = nt_demean(x);
x = (x .* nn) + mn;

x  = nt_fold(x, nsample);
w  = nt_fold(w, nsample);
ww = nt_fold(ww, nsample);

x(:, iNan, :)  = nan;
x(:, iZero, :) = 0;

if VERBOSE > 0, disp(['Final power ratio: ', num2str(nt_wpwr(x) / p0)]); end
end