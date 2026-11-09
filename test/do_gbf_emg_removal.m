function EEG = do_gbf_emg_removal(EEG, cfg)
% DO_GBF_EMG_REMOVAL Suppresses broadband EMG noise in continuous EEG data
% by projecting high-frequency activity onto a smooth cortical subspace.
%
% Inputs:
%   EEG       : EEGLAB data structure
%   K_inside  : Leadfield matrix for valid brain sources (channels x sources)
%   leadfield : Structure containing .pos, .tri, and .inside fields
%   cfg       : Configuration structure with fields:
%               .freq_cutoff (default: 20 Hz)
%               .num_modes   (default: 50 modes per hemisphere)
%               .beta        (default: 1e-4)
%
% Outputs:
%   EEG       : EEGLAB data structure with cleaned .data

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'freq_cutoff'), cfg.freq_cutoff = 20; end
if ~isfield(cfg, 'num_modes'),   cfg.num_modes   = 30; end
if ~isfield(cfg, 'beta'),        cfg.beta        = 0.05; end

fprintf('\n================================\n');
fprintf('GBF EMG Subspace Projection\n');
fprintf('================================\n');

% 0. Leadfield Preparation & Covariance Calculation
Table1 = readtable('C:\DATA\MATLAB\EEG\2_OTHER_DATA\FULL_CLINICAL_TABLE_2026-10-13.txt', 'VariableNamingRule', 'preserve');
Table1.ID = categorical(Table1.ID);
age = Table1.AGE(Table1.ID ==  EEG.ALSUTRECHT.subject.id);
age = mean(age);
fprintf('Average age: %.1f\n', age);

age_lf = 20:10:80;
[a, b] = min(abs(age_lf - age));
best_lf = age_lf(b);
% best_lf = 20;

file_name = ['leadfield_vrtx_' num2str(best_lf) 'yo.mat'];
file_path = fullfile('C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\template\eeg\BEM\age', file_name);
load(file_path, "leadfield");

leadfield = orient_leadfield(leadfield);
K_inside = cell2mat(leadfield.leadfield(leadfield.inside).');

indx = find(leadfield.inside, 1, "first");
num_ori = size(leadfield.leadfield{indx}, 2);
fprintf('Leadfield: %s\n', file_name);
fprintf('Using %d/%d dipoles from the leadfield.\n', size(K_inside, 2), num_ori*sum(leadfield.inside));

% 1. Extract inverse operator using the GBF MAP function
fprintf('1. Calculating GBF spatial filter (%d modes, beta=%g)...\n', cfg.num_modes, cfg.beta);
[W_inv, ~, ~] = do_gbf_map(K_inside, leadfield, cfg.num_modes, cfg.beta);

% 2. Reconstruct the spatial projection matrix (P_brain)
fprintf('2. Constructing sensor-to-sensor projection matrix...\n');
row_norms = vecnorm(K_inside, 2, 2);
K_norm = K_inside ./ row_norms;

% Apply diagonal scaling to correct for the row normalisation done inside do_gbf_map
P_brain = diag(row_norms) * K_norm * W_inv * diag(1 ./ row_norms);

% 3. Frequency split
fprintf('3. Splitting data at %d Hz...\n', cfg.freq_cutoff);
X = double(EEG.data);

% MATLAB's highpass operates along the first dimension (time), so we transpose
X_HF = highpass(X', cfg.freq_cutoff, EEG.srate)';

% Exact subtraction preserves the phase of the low-frequency data perfectly
X_LF = X - X_HF;

% 4. Apply projection to high-frequency data
fprintf('4. Applying cortical projection to high-frequency band...\n');
X_HF_clean = P_brain * X_HF;

% 5. Recombine
fprintf('5. Recombining low and cleaned high frequencies...\n');
% Cast back to single precision to adhere to standard EEGLAB formatting
data_clean = X_LF + X_HF_clean;

% Check
EEG_NEW = EEG;
EEG_NEW.data = data_clean;
cut_off = 1;
EEG_PLOT     = do_filtering_method(EEG, cut_off, cfg);
EEG_NEW_PLOT = do_filtering_method(EEG_NEW, cut_off, cfg);
vis_artifacts(EEG_NEW_PLOT, EEG_PLOT);

% Store
EEG.data = data_clean;

% Append a log note to the EEGLAB history
update_str = sprintf('GBF EMG removal applied: >%d Hz, %d modes, beta=%g', cfg.freq_cutoff, cfg.num_modes, cfg.beta);
if isfield(EEG, 'comments') && ischar(EEG.comments)
    EEG.comments = char(EEG.comments, update_str);
else
    EEG.comments = update_str;
end

fprintf('GBF cleaning complete.\n');
end

function [W, A_inside, lambda] = do_gbf_map(K_inside, leadfield, num_modes, beta)
% GBF_MAP Calculates a geometry-aware spatial filter using the
% Geometric Basis Function (GBF) framework and MAP estimation.
% Code:
% https://github.com/ncclab-sustech/GBFs/blob/master/utils/seeg_with_eeg.py
% https://github.com/ncclab-sustech/GBFs/blob/master/demo_cognitive.py

% 1. Geometry Extraction & Unit Conversion
% K_inside is in Ohm/mm
assert(strcmpi(leadfield.unit, 'mm'));
vertices = leadfield.pos; % mm
faces    = leadfield.tri;
inside   = leadfield.inside;

% 2. Row-wise Leadfield Normalisation (As seen in Python logic)
% Normalises each sensor (row) to a vector norm of 1.
row_norms = vecnorm(K_inside, 2, 2);
K_norm = K_inside ./ row_norms;

% 3. Split by Index (First half vs Second half)
n_total = size(vertices, 1);
n_half = n_total / 2;

if mod(n_total, 2) ~= 0
    error('Vertex count is not even; cannot split into equal halves.');
end

idx1 = 1:n_half;
idx2 = (n_half + 1):n_total;

% Assign names based on the mean X-coordinate
mean_x1 = mean(vertices(idx1, 1));
mean_x2 = mean(vertices(idx2, 1));

if mean_x1 < mean_x2
    fprintf('Detected L/R indices: 1:%d are left, %d:%d are right.\n', n_half, n_half+1, n_total);
else
    fprintf('Detected L/R indices: 1:%d are right, %d:%d are left.\n', n_half, n_half+1, n_total);
end

% % Check
% check_hemi = zeros(size(vertices, 1), 1);
% check_hemi(idx1) = -1; % LH will be one colour
% check_hemi(idx2) = 1;  % RH will be another colour
% figure;
% ft_plot_mesh(leadfield, 'vertexcolor', check_hemi);
% view([-90 0]); % View from the left side
% title('Hemisphere Split: Blue = LH?, Yellow = RH?');
% colorbar;

% 4. Extract Faces for each half
% Find faces where all three vertices belong to the same half
faces1_mask = all(ismember(faces, idx1), 2);
faces2_mask = all(ismember(faces, idx2), 2);

faces1 = faces(faces1_mask, :);
faces2 = faces(faces2_mask, :) - n_half; % Re-index for local Laplacian

% 5. Laplacian & Basis Extraction (num_modes per hemisphere)
fprintf('Processing left hemisphere...\n');
[S1, M1] = construct_cotangent_laplacian(vertices(idx1, :), faces1);

fprintf('Processing right hemisphere...\n');
[S2, M2] = construct_cotangent_laplacian(vertices(idx2, :), faces2);

[A1, D1] = eigs(S1, M1, num_modes, 'sm');
[A2, D2] = eigs(S2, M2, num_modes, 'sm');

% 6. Eigenvalue Extraction & Singularity Fix
lambda1 = abs(diag(D1));
lambda2 = abs(diag(D2));

% Handle the first trivial zero-mode by borrowing Mode 2's eigenvalue
lambda1(1) = lambda1(2);
lambda2(1) = lambda2(2);

% Sort eigenvalues in ascending order just in case eigs returned them unordered
[lambda1, sort_idx1] = sort(lambda1, 'ascend');
A1 = A1(:, sort_idx1);

[lambda2, sort_idx2] = sort(lambda2, 'ascend');
A2 = A2(:, sort_idx2);

% NORMALISE TO (0, 1):
% Scale so the highest mode sits at ~0.95 to keep -log(lambda) stable and positive
lambda1 = lambda1 / (max(lambda1) * 1.05);
lambda2 = lambda2 / (max(lambda2) * 1.05);
lambda  = [lambda1; lambda2];

assert(all(lambda < 1 & lambda > 0), 'Eigenvalues must be strictly between 0 and 1');
fprintf('Normalised eigenvalue spectrum: min = %.4e, max = %.4e\n', min(lambda), max(lambda));

% 7. Block-Diagonal Basis Construction
% This keeps the hemispheres separate in the solution space
A_full = zeros(n_total, num_modes * 2);
A_full(idx1, 1:num_modes) = A1;
A_full(idx2, num_modes + 1:end) = A2;

% % Check a Left Hemisphere Mode (should be between 1 and num_modes)
% mode_idx = 3;
% figure;
% subplot(1,2,1);
% ft_plot_mesh(leadfield, 'vertexcolor', A_full(:, mode_idx));
% title(['Basis Function ', num2str(mode_idx), ' (Should be LH Only)']);
% % Check a Right Hemisphere Mode (should be between num_modes+1 and 2*num_modes)
% subplot(1,2,2);
% ft_plot_mesh(leadfield, 'vertexcolor', A_full(:, num_modes + mode_idx));
% title(['Basis Function ', num2str(num_modes + 2), ' (Should be RH Only)']);

% Subset the Basis matrix to the "inside" sources
A_inside = A_full(inside, :);

% 8. MAP Estimation
% Effective leadfield in basis space
L = K_norm * A_inside;
L_sq = L' * L;

% Logarithmic Prior (Matches Python logic: -1/log(lambda))
Sigma_inv = diag(-beta ./ log(lambda));

fprintf('Computing the MAP estimator...\n');
fprintf('Ratio L_sq to Prior: %.2f\n', trace(L_sq) / trace(Sigma_inv)); % Around 5? Or 1-5?
W_theta = (L_sq + Sigma_inv) \ L';

% 9. Output Construction
% Final filter W maps sensors to source vertices.
W = A_inside * W_theta;

end

function [S, M] = construct_cotangent_laplacian(vertices, faces)
% CONSTRUCT_COTANGENT_LAPLACIAN Creates the stiffness and mass matrices
% for a three-dimensional triangular mesh.

num_vertices = size(vertices, 1);
num_faces = size(faces, 1);

% Initialise arrays for sparse matrix assembly.
% We allocate 6 entries per face (2 for each of the 3 edges).
I = zeros(num_faces * 6, 1);
J = zeros(num_faces * 6, 1);
V = zeros(num_faces * 6, 1);

% 1. Calculate the edge vectors for every triangle in the mesh.
% v1 goes from vertex 1 to vertex 2.
v1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
% v2 goes from vertex 2 to vertex 3.
v2 = vertices(faces(:, 3), :) - vertices(faces(:, 2), :);
% v3 goes from vertex 3 to vertex 1.
v3 = vertices(faces(:, 1), :) - vertices(faces(:, 3), :);

% 2. Calculate the area of each triangle.
% The magnitude of the cross product of two edge vectors gives twice
% the area of the triangle.
cross_prod = cross(v1, -v3, 2);
triangle_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));

% 3. Calculate the cotangents of the internal angles.
% We use the geometric identity: cot(theta) = dot(u,v) / norm(cross(u,v)).
cot1 = dot(-v3, v1, 2) ./ sqrt(sum(cross(-v3, v1, 2).^2, 2));
cot2 = dot(-v1, v2, 2) ./ sqrt(sum(cross(-v1, v2, 2).^2, 2));
cot3 = dot(-v2, v3, 2) ./ sqrt(sum(cross(-v2, v3, 2).^2, 2));

% 4. Assemble the stiffness matrix S (the cotangent weights).
idx = 1;
for i = 1:3
    i1 = i;
    i2 = mod(i, 3) + 1;
    i3 = mod(i + 1, 3) + 1;

    edges1 = faces(:, i2);
    edges2 = faces(:, i3);

    if i == 1
        cot_weight = cot1;
    elseif i == 2
        cot_weight = cot2;
    else
        cot_weight = cot3;
    end

    % Add the symmetric off-diagonal weights for the current edge.
    I(idx:idx+num_faces-1) = edges1;
    J(idx:idx+num_faces-1) = edges2;
    V(idx:idx+num_faces-1) = -0.5 * cot_weight;
    idx = idx + num_faces;

    I(idx:idx+num_faces-1) = edges2;
    J(idx:idx+num_faces-1) = edges1;
    V(idx:idx+num_faces-1) = -0.5 * cot_weight;
    idx = idx + num_faces;
end

% Construct the sparse matrix from the triplets.
S = sparse(I(1:idx-1), J(1:idx-1), V(1:idx-1), num_vertices, num_vertices);

% The diagonal of the stiffness matrix must equal the negative sum
% of the off-diagonal elements in the corresponding row.
S = S - spdiags(sum(S, 2), 0, num_vertices, num_vertices);

% 5. Assemble the lumped mass matrix M.
% This represents the local area associated with each vertex.
% We distribute one third of each triangle's area to its three vertices.
vertex_areas = accumarray(faces(:), repmat(triangle_areas / 3, 3, 1), [num_vertices, 1]);
M = spdiags(vertex_areas, 0, num_vertices, num_vertices);
end

function lf_out = orient_leadfield(lf_in)
% ORIENT_LEADFIELD Projects unconstrained leadfields onto surface normals.
%
% Input:
%   lf_in  - FieldTrip leadfield structure with fields:
%            .leadfield {Npos x 1 cell} of [Nchan x 3]
%            .ori       {Npos x 1 cell} of [1 x 3] or [3 x 1] (unit surface normals)
%            .inside    [Npos x 1 logical]
%
% Output:
%   lf_out - Leadfield structure with oriented leadfields [Nchan x 1]

lf_out = lf_in;

inside_idx = find(lf_in.inside);

for i = 1:numel(inside_idx)
    pos_idx = inside_idx(i);

    L = lf_in.leadfield{pos_idx}; % [Nchan x 3]
    u = lf_in.ori{pos_idx};       % [1 x 3] or [3 x 1]

    if isempty(L) || isempty(u)
        continue;
    end

    % Ensure column vector [3 x 1] and unit length
    u = u(:);
    norm_u = norm(u);
    if norm_u > 0
        u = u / norm_u;
    end

    % Project: [Nchan x 3] * [3 x 1] -> [Nchan x 1]
    lf_out.leadfield{pos_idx} = L * u;
end

% Update dimord metadata
if isfield(lf_out, 'leadfielddimord')
    lf_out.leadfielddimord = '{pos}_chan';
end
end

% function plot_cortical_mode(vertices, faces, mode_data, mode_number)
% % PLOT_CORTICAL_MODE Visualises a specific geometric basis function
% % on the 3D cortical surface mesh for debugging.
% %
% % Inputs:
% %   vertices    : Spatial coordinates of the cortical surface (Sources by 3)
% %   faces       : Triangulation connectivity of the surface (Triangles by 3)
% %   mode_data   : The spatial pattern of the mode to plot (Sources by modes)
% %   mode_number : Integer identifying the mode number for the title
%
% mode_data = mode_data(:, mode_number);
%
% % Create a new figure with a white background
% figure('Name', sprintf('Cortical Mode %d', mode_number), 'Color', 'w');
%
% % Plot the 3D mesh using trisurf
% % 'FaceColor', 'interp' smoothly interpolates the colours across the triangles
% trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
%     mode_data, 'EdgeColor', 'none', 'FaceColor', 'interp');
%
% % Optimise the camera view and axis proportions
% axis equal;
% axis off;
% view([-90, 0]); % Default to a lateral view; adjust angles as needed
%
% % Add lighting and shading to enhance the 3D curvature of the brain folds
% camlight('headlight');
% lighting gouraud;
% material dull;
%
% % Add a colourbar and set the colour map
% % 'jet' or standard diverging colour maps work well for oscillating modes
% colorbar;
% colormap(gca, 'jet');
%
% % Centre the colour limits around zero
% % This ensures that negative values and positive values are symmetrical
% max_val = max(abs(mode_data));
% if max_val > 0
%     clim([-max_val, max_val]); % Use caxis([-max_val, max_val]) for older MATLAB versions
% end
%
% % Add an informative title
% title(sprintf('Geometric Basis Function: Mode %d', mode_number), 'FontSize', 14);
% end

% function areas = calculate_face_areas(vertices, faces)
% % CALCULATE_FACE_AREAS Computes the surface area of each triangular face.
% %
% % Inputs:
% %   vertices : [N x 3] matrix of vertex coordinates (in metres)
% %   faces    : [M x 3] matrix of vertex indices for each triangle
% %
% % Outputs:
% %   areas    : [M x 1] vector containing the area of each face
%
% % 1. Extract the coordinates for the three vertices of every triangle
% A = vertices(faces(:, 1), :);
% B = vertices(faces(:, 2), :);
% C = vertices(faces(:, 3), :);
%
% % 2. Compute two edge vectors for each triangle
% vec1 = B - A;
% vec2 = C - A;
%
% % 3. The magnitude of the cross product of two vectors
% % equals the area of the parallelogram they form.
% cp = cross(vec1, vec2, 2);
%
% % 4. Half the magnitude of the cross product is the triangle area.
% % This is vectorised for speed.
% areas = 0.5 * sqrt(sum(cp.^2, 2));
% end