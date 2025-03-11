function [M, y, score, prop] = nt_sca(x,ncomp)
% [M,y,score,proportion] = nt_sca(x,ncomp) - shared component analysis
%
% x: data (time X channels) (can be cell array)
% ncomp: keep only ncomp components (faster)
%
% M: SCA transform matrix
% z: transformed data
% score: sharedness score, per component
% prop: proportion of power accounted for, per component

% THRESH=10^-12;

%{
Todo: 
- allow xvalidation
- operate on covariance matrices
%}

% SDukic, edited February 2025

if nargin<1; error('!'); end
if nargin<2; ncomp=[]; end
if isempty(ncomp); ncomp = size(x,2); end

if diff(size(x)) ~= 0
    isRaw = true;
    if iscell(x)
        % xx = [];
        % for iTrial = 1:numel(x)
        %     xx = [xx; x{iTrial}];
        % end
        % x = xx; clearvars xx;
        isCell = true;
        [NPTS, NCHN] = size(x{1}); NTRL = length(x);
        x = cat(1,x{:});
    else
        isCell = false;
    end
    x = nt_demean(x);
    C0 = nt_cov(x);  % initial covariance
else
    isRaw = false;
    isCell = false;
    C0 = x;
end

T = eye(size(x,2)); % current transform
M = eye(size(x,2)); % result

score = NaN(ncomp,1);
for iComp = 1:ncomp
    C = T' * C0 * T;                 % current covariance
    N = diag(1./sqrt(diag(C)));      % normalizing matrix
    N(isnan(N)) = 0;
    C = N' * C * N;                  % normalize current covariance
    [topcs, ev] = nt_pcarot(C);      % PCA
    frompcs = pinv(topcs);
    M(:,iComp) = T * N * topcs(:,1); % keep first PC
    T = T * N * (topcs(:,2:end) * frompcs(2:end,:)); % project out first PC
    score(iComp) = ev(1);
end

if ncomp < size(x,2)
    % fill rest of transform matrix with leftover (unprocessed) dimensions
    C = T' * C0 * T;            % current covariance
    N = diag(1./sqrt(diag(C))); % normalizing matrix
    N(isnan(N)) = 0;
    C = N' * C * N;             % normalize current covariance
    topcs = nt_pcarot(C);       % PCA
    T = T * topcs;
    M(:,ncomp+1:end) = T(:,1:(size(x,2)-ncomp));
end

prop = diag(M'*C*M);

% figure(1); clf; plot(prop); pause;
if nargout > 1
    if isRaw
        y = x * M;
    else
        y = M' * x * M;
    end
end

if isCell
    y = reshape(y',NCHN,NPTS,NTRL);
    yy = cell(1,NTRL);
    for iTrial = 1:NTRL
        yy{iTrial} = y(:,:,iTrial)';
    end
    y = yy;
end

% % test code
% if 0
%     x=randn(1000,10);
%     [M,y]=nt_sca(x);
% end
%
% if 0
%     % data are 11 chan:
%     % 10 chan share same source (sine),
%     % 1 chan is different source (noise) with higher variance
%     x=randn(1000,10);
%     s=sin(2*pi*(1:1000)'/1000);
%     x=bsxfun(@plus,x,s); % add same to all
%     x=[x,10*randn(1000,1)]; % extra channel with large variance
%     %[y,M,score]=nt_sca_old(x);
%     MM=nt_sca(x); y=x*MM;
%     yy=nt_pca(x);
%     figure(1); clf;  plot(nt_normcol(s)'*nt_normcol(y)/size(s,1));
%     hold on; plot(nt_normcol(s)'*nt_normcol(yy)/size(s,1));    legend('sca','pca');
% end
%
%
% if 0
%     % two shared sources
%     x=randn(1000,10);
%     s=sin(2*pi*(1:1000)'/1000);
%     s2=sin(2*pi*2*(1:1000)'/1000);
%     x=x+s*rand(1,10); % add same to all
%     x=x+s2*rand(1,10); % add same to all
%     x=[x,10*randn(1000,3)]; % extra channel with large variance
%     %[y,M,score]=nt_sca_old(x);
%     MM=nt_sca(x); yyy=x*MM;
%     yy=nt_pca(x);
%     figure(1); clf;
%     subplot 121; %bar(abs(nt_normcol(s)'*nt_normcol(y)/size(s,1)));
%     hold on; bar(abs(nt_normcol(s)'*nt_normcol(yyy)/size(s,1)));
%     hold on; bar(abs(nt_normcol(s)'*nt_normcol(yy)/size(s,1)));    legend('sca','pca'); title('source 1');
%     subplot 122; %bar(abs(nt_normcol(s2)'*nt_normcol(y)/size(s,1)));
%     hold on;bar(abs(nt_normcol(s2)'*nt_normcol(yyy)/size(s,1)));
%     hold on; bar(abs(nt_normcol(s2)'*nt_normcol(yy)/size(s,1)));    legend('sca','pca'); title('source 2');
% end