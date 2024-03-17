function rank2 = get_rank(data)

% EEGLAB sometimes forces data to be single
data = double(data);

% Reshape if data is cut into trials
if size(data,3)>1
    data = reshape(data,size(data,1),[]);
end

% MATLAB's func
rank1 = rank(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here: alternate computation of the rank by Sven Hoffman
% tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
% [E, D] = eig(cov(double(tmpdata')));
% rank2 = sum(diag(D) > 1e-6);

% This threshold might not work for MEG data???
if ~diff(size(data))
    % if  sqaure matrix
    rank2 = sum(eig(data) > 1e-7);
else
    rank2 = sum(eig(cov(data')) > 1e-7);
end

if rank1 ~= rank2
    % fprintf('Warning: fixing rank computation inconsistency (%d vs %d)
    % most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
    % rank2 = max(rank1, rank2);

    fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because of very correlated channels...\n', rank1, rank2);
    rank2 = min(rank1, rank2);
end

end