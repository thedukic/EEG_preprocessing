function data = spatial_filtering(data,elecpos,neighbours,Nneighbours)
% Nnb - number of neighbours, if 0, then it uses all neighbours
%
% Baseline correction is probably good to do beforehand
% This method is problematic if there are multiple channels affected at once
%
% Add:
% - Check if each channel has at least Nnb
%

[NCHN, NTRL] = size(data);

for i = 1:NCHN
    % Neighbour list
    ns = find(logical(neighbours(i,:)));

    % Remove the current electrode from the neighbour list
    ns(ns==i) = [];

    % ds = sqrt(sum((elecpos(ns,:) - repmat(elecpos(i,:), numel(ns), 1)).^2, 2));
    ds = vecnorm((elecpos(ns,:)-elecpos(i,:))')';
    assert(length(ns)==length(ds));

    if Nneighbours>0
        ind = 1:length(ns)-Nneighbours;
        [~, tmp] = sort(ds,'descend');
        ns(tmp(ind)) = [];
        ds(tmp(ind)) = [];
        assert(length(ns)==Nneighbours);
    end

    [x, indx] = sort(data(ns,:));
    ws = repmat([1; (1./ds)],1,size(x,2));
    ws = ws(indx);

    x([1 end],:) = [];
    ws([1 end],:) = [];

    ws = ws./sum(ws);
    data(i,:) = sum(x.*ws);

    % Smooth over time
    % data(i,:) = smoothdata(data(i,:),"sgolay",10);
end

end