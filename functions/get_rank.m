function rank2 = get_rank(data)

% EEGLAB may force data to be single
data = double(data);

% Reshape if data is cut into trials
if size(data,3)>1
    data = reshape(data,size(data,1),[]);
end

% MATLAB function
rank1 = rank(data);

% Alternative estiamte of the rank by Sven Hoffman
% Avoids false high ranks although some electrodes were interpolated
if ~diff(size(data))
    % If it is a square/covariance matrix
    rank2 = sum(eig(data) > 1e-7);
else
    rank2 = sum(eig(cov(data.')) > 1e-7);
end

% Check
if rank1 ~= rank2
    fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because of very correlated channels...\n', rank1, rank2);
    rank2 = min(rank1, rank2);
end

end