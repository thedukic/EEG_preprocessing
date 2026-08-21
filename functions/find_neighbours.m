function neighbours = find_neighbours(data, target_region)
% FIND_NEIGHBOURS_NEW Calculates spatial or topological neighbours for TFCE.
%
% Usage:
%   neighbours = find_neighbours(locs_struct)
%   neighbours = find_neighbours(source_struct, region_to_plot)
%
% Inputs:
%   data  - For 'sensor': Structure array with .X, .Y, .Z fields (eLoc).
%           For 'source': Structure with three fields:
%               .tri    (F x 3 matrix of cortical mesh triangles)
%               .tissue (V x 1 numeric vector mapping each vertex to an atlas region ID)
%               .pos    (V x 3 matrix of vertex coordinates, required for plotting)
%   target_region - (Optional) Integer specifying which source region to highlight
%                   in the plot. Defaults to 1.

if nargin < 2
    target_region = 10;
end

do_plot = false;

if length(data) > 1
    level = 'sensor';
else
    level = 'source';
end

if strcmpi(level, 'sensor')
    % Extract coordinates
    x = [data(:).X]';
    y = [data(:).Y]';
    z = [data(:).Z]';
    vertices = [x, y, z];
    nRegions = size(vertices, 1);

    % 1. Flat projection for border detection
    z2 = z - max(z);
    hypotxy = hypot(x, y);
    R = hypot(hypotxy, z2);
    PHI = atan2(z2, hypotxy);
    TH = atan2(y, x);
    PHI(PHI < 0.001) = 0.001;
    R2 = R ./ cos(PHI) .^ .2;
    X = R2 .* cos(TH);
    Y = R2 .* sin(TH);

    % 2. Optimize Head Center using fminsearch
    mass = mean(vertices, 1);
    diffvert = bsxfun(@minus, vertices, mass);
    R0 = mean(sqrt(sum(diffvert.^2, 2)));

    % Anonymous function replicating dist_sph
    objFun = @(v) mean(abs(sqrt(sum(bsxfun(@minus, vertices, v(1:3)).^2, 2)) - v(4)));
    vec0 = [mass, R0];
    minn = fminsearch(objFun, vec0);
    HeadCenter = minn(1:3);

    % 3. Project to sphere using optimized center
    coordC = bsxfun(@minus, vertices, HeadCenter);
    coordC = bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2, 2)));

    % Tesselation of the sensor array
    faces = convhulln(coordC);

    % 4. Remove artificial triangles capping the bottom of the net
    border = convhull(X, Y);
    iInside = ~(ismember(faces(:,1), border) & ismember(faces(:,2), border) & ismember(faces(:,3), border));
    faces = faces(iInside, :);

    % 5. Threshold based on perimeter
    v1 = vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 2), :);
    v3 = vertices(faces(:, 3), :);
    triPerimeter = sqrt(sum((v1 - v2).^2, 2)) + ...
        sqrt(sum((v2 - v3).^2, 2)) + ...
        sqrt(sum((v3 - v1).^2, 2));
    thresholdPerim = mean(triPerimeter) + 3 * std(triPerimeter);
    faces(triPerimeter > thresholdPerim, :) = [];

    % 6. Build neighbour lists
    output = cell(nRegions, 1);
    for n = 1:nRegions
        [r, ~] = find(faces == n);
        output{n} = unique(faces(r, :));
        output{n} = reshape(output{n}, 1, []); % Ensure row vector
    end

elseif strcmpi(level, 'source')
    % Extract atlas mesh and labels
    faces = data.tri;
    labels = data.tissue;

    % Ensure labels is a column vector
    labels = labels(:);

    % Extract all unique edges from the mesh (vertex pairs)
    edges = [faces(:, 1), faces(:, 2); ...
        faces(:, 2), faces(:, 3); ...
        faces(:, 3), faces(:, 1)];

    % Map vertex edges to atlas region edges
    region_edges = [labels(edges(:, 1)), labels(edges(:, 2))];

    % Filter out connections within the same region and unlabelled vertices (e.g., 0 or NaN)
    valid_idx = (region_edges(:, 1) ~= region_edges(:, 2)) & ...
        (region_edges(:, 1) > 0) & (region_edges(:, 2) > 0) & ...
        ~isnan(region_edges(:, 1)) & ~isnan(region_edges(:, 2));
    valid_borders = region_edges(valid_idx, :);

    % Make bidirectional and find unique pairs
    all_pairs = unique([valid_borders; valid_borders(:, 2), valid_borders(:, 1)], 'rows');

    % Number of regions based on the maximum ID found in the atlas
    nRegions = max(labels(~isnan(labels) & labels > 0));
    output = cell(nRegions, 1);
    for n = 1:nRegions
        % Find regions connected to region 'n'
        region_neighbours = all_pairs(all_pairs(:, 1) == n, 2);
        % Include the region itself in the list
        output{n} = unique([n; region_neighbours])';
    end
else
    error('Level must be specified as either ''sensor'' or ''source''.');
end

% Convert cell array to padded matrix format required by TFCE
nz = max(cellfun(@numel, output));
neighbours = zeros(nRegions, nz);
for i = 1:nRegions
    neighbours(i, 1:length(output{i})) = output{i};
end

% --- PLOTTING SECTION ---
if do_plot
    figure('Color','w','Position',[50, 50, 600, 500]);
    axes('Color','w');

    if strcmpi(level, 'sensor')
        patch('Vertices', vertices, 'Faces', faces, ...
            'FaceVertexCData', repmat([1 1 1], [length(vertices), 1]), ...
            'Marker', 'o', 'LineWidth', 1, ...
            'FaceColor', [.5 .5 .5], 'FaceAlpha', .9, ...
            'EdgeColor', [0 0 0], 'EdgeAlpha', 1, ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', 'flat', ...
            'MarkerSize', 12, 'BackfaceLighting', 'lit');
        title('Sensor Triangulation');

    elseif strcmpi(level, 'source')
        % To plot the source, coordinate data must be provided
        if isfield(data, 'pos')
            src_vertices = data.pos;
        elseif isfield(data, 'vertices')
            src_vertices = data.vertices;
        else
            warning('Cannot plot source model. Please include .pos or .vertices in your data structure.');
            return;
        end

        % Initialise all vertices to grey
        vertex_colors = repmat([0.75 0.75 0.75], size(src_vertices, 1), 1);

        if target_region > 0 && target_region <= nRegions
            % Extract the neighbours for the target region, excluding the target itself
            neigh_regions = setdiff(output{target_region}, target_region);
            num_neigh = length(neigh_regions);

            % Generate distinct colours for each neighbour using the 'lines' colormap
            neigh_cmap = lines(num_neigh);

            % Assign each neighbour its unique colour
            for k = 1:num_neigh
                idx_neigh = (labels == neigh_regions(k));
                vertex_colors(idx_neigh, :) = repmat(neigh_cmap(k, :), sum(idx_neigh), 1);
            end

            % Set colour for the target region (Red)
            idx_target = (labels == target_region);
            vertex_colors(idx_target, :) = repmat([0.9 0.2 0.2], sum(idx_target), 1);

            title(sprintf('Source Region %d (Red) and its %d Neighbours (Coloured)', target_region, num_neigh));
        else
            title('Source Atlas Regions (Target Region Out of Bounds)');
        end

        % Plot the brain mesh using the custom RGB vertex colours
        patch('Vertices', src_vertices, 'Faces', faces, ...
            'FaceVertexCData', vertex_colors, ...
            'FaceColor', 'interp', ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1, ...
            'BackfaceLighting', 'lit');
    end

    % Common lighting and view settings
    material([0.5 0.50 0.20 1.00 0.5]);
    camlight;
    lighting phong;
    view(48,15);
    axis equal;
    axis off;
    cameratoolbar('Show');
    rotate3d on;
end

end