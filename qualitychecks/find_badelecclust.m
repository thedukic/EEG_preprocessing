function [clusters, cluster_sizes, distance_threshold] = find_badelecclust(bad_electrode_indices, all_electrode_coords, distance_threshold)
%FIND_SPATIAL_ELECTRODE_CLUSTERS Identifies clusters based on 3D spatial proximity.
%
%   This function groups bad electrodes into clusters if they are within a
%   specified Euclidean distance of each other in XYZ space.
%
%   [clusters, cluster_sizes, distance_threshold] = FIND_SPATIAL_ELECTRODE_CLUSTERS(bad_electrode_indices, all_electrode_coords, distance_threshold)
%
%   INPUTS:
%   - bad_electrode_indices: A vector of indices for the bad electrodes.
%   - all_electrode_coords:  An N-by-3 matrix where N is the total number
%                            of electrodes, and each row contains the [X, Y, Z]
%                            coordinates for an electrode.
%   - distance_threshold:    (Optional) A scalar value. Two electrodes are
%                            considered neighbours if the distance between them
%                            is less than or equal to this threshold.
%                            If this argument is omitted or passed as an empty
%                            array `[]`, the threshold is calculated automatically.
%
%   OUTPUTS:
%   - clusters:      A cell array where each cell contains a vector of
%                    electrode indices belonging to one spatial cluster.
%   - cluster_sizes: A vector containing the size of each corresponding
%                    cluster in the 'clusters' cell array.
%   - distance_threshold: The threshold value used for clustering (either
%                         the value passed in or the one automatically calculated).
%
%   Example Usage:
%       all_coords = rand(64, 3) * 100; % 64 electrodes
%       bad_indices = [1, 5, 6, 23, 24, 25, 40];
%
%       % Let the function determine the threshold by omitting the argument
%       [found_clusters, sizes, used_thresh] = find_spatial_electrode_clusters(bad_indices, all_coords);


% --- Input Validation ---
if isempty(bad_electrode_indices)
    clusters = {};
    cluster_sizes = [];
    % Set threshold to NaN if no bad electrodes are provided.
    if nargin < 3 || isempty(distance_threshold)
        distance_threshold = NaN;
    end
    return;
end

if size(all_electrode_coords, 2) ~= 3
    error('The `all_electrode_coords` matrix must have 3 columns for X, Y, and Z.');
end

% --- Threshold Calculation ---
% Check if the distance_threshold argument was not provided or was passed as empty.
if nargin < 3 || isempty(distance_threshold)
    % Automatically determine the threshold based on electrode geometry.

    % Calculate pairwise distances between ALL electrodes.
    all_dist_matrix = pdist2(all_electrode_coords, all_electrode_coords);

    % To find the nearest neighbour for each electrode, we ignore the
    % distance to itself (which is 0) by setting the diagonal to infinity.
    all_dist_matrix(logical(eye(size(all_dist_matrix)))) = inf;

    % Find the minimum distance for each electrode to its nearest neighbour.
    min_distances = min(all_dist_matrix, [], 2);

    % Sort these minimum distances to find the smallest ones.
    sorted_min_distances = sort(min_distances);

    % Determine the number of electrodes that make up the bottom 20%.
    num_to_average = ceil(0.15 * numel(sorted_min_distances));

    % Calculate the threshold as the average of these smallest distances.
    distance_threshold = mean(sorted_min_distances(1:num_to_average));
    distance_threshold = 1.5 * distance_threshold;

    % fprintf('Automatically determined distance threshold: %.2f\n', distance_threshold);
end

% --- Main Logic ---

% Get the coordinates of only the bad electrodes.
bad_coords = all_electrode_coords(bad_electrode_indices, :);
num_bad_electrodes = numel(bad_electrode_indices);

% Calculate the pairwise Euclidean distance between all bad electrodes.
% This creates a square matrix where dist_matrix(i, j) is the
% distance between bad electrode i and bad electrode j.
dist_matrix = pdist2(bad_coords, bad_coords);

% Create an adjacency matrix. Two electrodes are "adjacent" (connected)
% if their distance is within the threshold. Here we include the diagonal
% (distance 0) to ensure single-electrode clusters are handled correctly
% by the graph traversal logic.
adjacency_matrix = (dist_matrix <= distance_threshold);

% --- Find Connected Components (Clusters) using Graph Traversal ---

% Keep track of electrodes that have already been assigned to a cluster.
visited = false(1, num_bad_electrodes);
clusters = {};

for i = 1:num_bad_electrodes
    if ~visited(i)
        % This electrode has not been visited, so it starts a new cluster.

        % Use a Breadth-First Search (BFS) to find all connected electrodes.
        queue = i;          % Start the search with the current electrode.
        visited(i) = true;

        current_cluster_indices = []; % Stores indices within the bad_indices vector

        head = 1;
        while head <= numel(queue)
            current_node = queue(head);
            head = head + 1;

            current_cluster_indices = [current_cluster_indices, current_node];

            % Find all neighbours of the current node from the adjacency matrix.
            neighbours = find(adjacency_matrix(current_node, :));

            % Add unvisited neighbours to the queue.
            for j = 1:numel(neighbours)
                neighbour_node = neighbours(j);
                if ~visited(neighbour_node)
                    visited(neighbour_node) = true;
                    queue = [queue, neighbour_node];
                end
            end
        end

        % Convert the cluster indices back to the original electrode numbers.
        clusters{end+1} = bad_electrode_indices(sort(current_cluster_indices));
    end
end

% Calculate the sizes of the found clusters.
cluster_sizes = cellfun(@numel, clusters);

end
