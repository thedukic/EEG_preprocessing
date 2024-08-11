function rank2 = get_rank(data)
%
% Input data (Nchannels x Ntimepoints) or (Nchannels x Ntimepoints x Ntrials)
% EEGLAB may force data to be single
%

data = double(data);

% Reshape if data is cut into trials
if size(data,3)>1
    data = reshape(data,size(data,1),[]);
end

% MATLAB function
rank1 = rank(data);

% Alternative estimate of the rank by Sven Hoffman
% Avoids false high ranks although some electrodes were interpolated
if ~diff(size(data))
    % If it is a square/covariance matrix
    rank2 = sum(eig(data) > 1e-7);
else
    rank2 = sum(eig(cov(data.')) > 1e-7);
end

% Check
if rank1 ~= rank2
    fprintf('Fixing rank computation inconsistency (%d vs %d) most likely because of very correlated channels or spherical channel interpolation...\n', rank1, rank2);
    rank2 = min(rank1, rank2);
end

end