function roiNeighbours = find_neighbours(roilocs,R)
% Define the radius R for neighbouring detection
% List of brain region centres (in 3D space: x, y, z)

% Extract
x = [roilocs(:).X]';
y = [roilocs(:).Y]';
z = [roilocs(:).Z]';
regionCentres = [x, y, z];

% Number of regions
numRegions = size(regionCentres, 1);

% Initialise a cell array to store neighbours for each region
neighbours = cell(numRegions, 1);

% Loop through each region and calculate distances to others
for i = 1:numRegions
    % Compute distances from the i-th region to all other regions
    distances = sqrt(sum((regionCentres - regionCentres(i, :)).^2, 2));

    % Find indices of regions within distance R (excluding itself)
    % neighbourIndices = find(distances < R & distances > 0);

    % Sort them by distance
    % -> the first is 0 mm distance
    [distances, indxSorted] = sort(distances);
    neighbourIndices = distances < R;
    neighbourIndices = indxSorted(neighbourIndices);

    % Store neighbour indices for the i-th region
    neighbours{i} = neighbourIndices;
end

% % Display results
% for i = 1:numRegions
%     fprintf('Region %d neighbours: %s\n', i, mat2str(neighbours{i}));
% end

% Convert for TFCE
maxLengths = max(cellfun(@length, neighbours));
roiNeighbours = zeros(numRegions, maxLengths);

for i = 1:numRegions
    roiNeighbours(i, 1:length(neighbours{i})) = neighbours{i};
end

% figure; hold on;
% K = 8;
% ft_plot_mesh(regionCentres);
% ft_plot_mesh(regionCentres(K,:),'vertexcolor','r','vertexsize',40);
% neighTmp = roiNeighbours(K,:);
% neighTmp(neighTmp == 0) = [];
% ft_plot_mesh(regionCentres(neighTmp,:),'vertexcolor','b','vertexsize',30);

% figure;
% for i = 1:numRegions
%     cla; mytopoplot(ones(1,128),roiNeighbours(i,:),[],nexttile(1));
%     pause;
% end

end