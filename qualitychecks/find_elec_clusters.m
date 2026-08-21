function [clusters, cluster_sizes] = find_elec_clusters(bad_electrodes)
%FIND_ELEC_CLUSTERS Identifies spatial clusters based on 3D triangulation.
%
%   This function groups bad electrodes into clusters by first calculating a
%   topological triangulation of the entire scalp cap (using find_neighbours),
%   and then finding connected components among the bad electrodes.
%
%   [clusters, cluster_sizes] = find_elec_clusters(bad_electrode_indices)
%
%   INPUTS:
%   - bad_electrode_indices: A vector of indices for the bad electrodes.
%
%   OUTPUTS:
%   - clusters:      A cell array where each cell contains a vector of
%                    electrode indices belonging to one spatial cluster.
%   - cluster_sizes: A vector containing the size of each corresponding
%                    cluster in the 'clusters' cell array.

if isempty(bad_electrodes)
    clusters = {};
    cluster_sizes = [];
    return;
end

% Find the channel coordiantes
load('biosemi128_eeglab.mat', 'chanlocs');

% Organise channel indices
channel_labels = {chanlocs.labels};
num_total_elec = length(channel_labels);
if iscell(bad_electrodes)
    bad_electrodes = find(ismember(channel_labels, bad_electrodes));
elseif islogical(bad_electrodes)
    bad_electrodes = find(bad_electrodes);
end

% 2. Calculate topological neighbours (calls the local function below)
neighbours_padded = find_neighbours(chanlocs);

% 3. Build a full logical N-by-N adjacency matrix
adjacency_matrix = false(num_total_elec, num_total_elec);
for i = 1:num_total_elec
    % Extract valid neighbour indices (ignoring the 0-padding)
    valid_neighs = neighbours_padded(i, neighbours_padded(i, :) > 0);
    adjacency_matrix(i, valid_neighs) = true;
end

% Ensure matrix symmetry (in case of directional triangulation artifacts)
adjacency_matrix = adjacency_matrix | adjacency_matrix';

% 4. Extract the sub-graph for only the bad electrodes
bad_adj = adjacency_matrix(bad_electrodes, bad_electrodes);

% Ensure the diagonal is 0 (no self-loops in the graph theory functions)
bad_adj(logical(eye(size(bad_adj)))) = false;

% 5. Create an undirected graph and find connected components instantly
G = graph(bad_adj);
bins = conncomp(G);

% Group the results by cluster
num_clusters = max(bins);
clusters = cell(1, num_clusters);
cluster_sizes = zeros(1, num_clusters);

for i = 1:num_clusters
    % Map the graph nodes back to the original electrode indices
    members = find(bins == i);
    clusters{i} = bad_electrodes(members);
    cluster_sizes(i) = length(members);
end

% Sort clusters by size (largest first)
[cluster_sizes, sort_idx] = sort(cluster_sizes, 'descend');
clusters = clusters(sort_idx);

end
