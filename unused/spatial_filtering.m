function dataOut = spatial_filtering(dataIn,elecpos,neighbours,NNeighboursSet)
% Nnb - number of neighbours, if 0, then it uses all neighbours
%
% Baseline correction is probably good to do beforehand
% This method is problematic if there are multiple channels affected at once
%
% Add:
% - Check if each channel has at least Nnb
%

[NCHN, NTRL] = size(dataIn);
dataOut = NaN(size(dataIn));
NNeighboursMax = NaN(NCHN,1);

NNeighboursSet = NNeighboursSet+2;

for i = 1:NCHN
    % Neighbour list
    ns = find(logical(neighbours(i,:)));

    % Remove the current electrode from the neighbour list
    ns(ns==i) = [];
    NNeighboursMax(i) = length(ns);
    disp(NNeighboursMax(i));

    % ds = sqrt(sum((elecpos(ns,:) - repmat(elecpos(i,:), numel(ns), 1)).^2, 2));
    ds = vecnorm((elecpos(ns,:) - elecpos(i,:))')';
    assert(length(ns) == length(ds));

    if NNeighboursSet>0
        ind = 1:length(ns)-NNeighboursSet;
        [~, tmp] = sort(ds,'descend');
        ns(tmp(ind)) = [];
        ds(tmp(ind)) = [];
        assert(length(ns)==NNeighboursSet);
    end

    [x, indx] = sort(dataIn(ns,:));
    ws = repmat([1; (1./ds)],1,size(x,2));
    ws = ws(indx);

    x([1 end],:) = [];
    ws([1 end],:) = [];

    ws = ws ./ sum(ws);
    dataOut(i,:) = sum(x.*ws);

    % Smooth over time
    % data(i,:) = smoothdata(data(i,:),"sgolay",10);
end

end