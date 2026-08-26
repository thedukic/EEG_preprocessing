function W_rotated = rotate_topotemplate(W_orig, chanlocs, angle_deg)
% ROTATE_CHANNEL_WEIGHTS Rotates a vector of channel weights and reassigns
% them to the original channel locations using inverse mapping and interpolation.
%
% This achieves the effect of rotating the ENTIRE SPATIAL PATTERN of weights
% by the specified angle while keeping the channel indices fixed.
%
% Inputs:
%   W_orig    - 128x1 vector of original channel weights (the data values).
%   X_orig    - 128x1 vector of original channel X-coordinates.
%   Y_orig    - 128x1 vector of original channel Y-coordinates.
%   angle_deg - Angle of rotation for the PATTERN, in degrees (e.g., 30).
%               Positive angle_deg means clockwise rotation of the pattern.
%
% Output:
%   W_rotated - 128x1 vector of new weights assigned to the original channel indices.
%
% Example:
%   % Assume X_coords, Y_coords, and W_weights are defined (e.g., load your data)
%   % W_rotated = rotate_channel_weights(W_weights, X_coords, Y_coords, 45);
%   % topoplot(W_rotated, X_coords, Y_coords, 'maplimits', [-1 1]);

fprintf('Starting pattern rotation by %.1f degrees (clockwise).\n', angle_deg);

indices = 1:128;
x = [chanlocs(indices).X]';
y = [chanlocs(indices).Y]';
z = [chanlocs(indices).Z]';

z2 = z - max(z);

hypotxy = hypot(x,y);
R = hypot(hypotxy,z2);
PHI = atan2(z2,hypotxy);
TH = atan2(y,x);

% Remove the too small values for PHI
PHI(PHI < 0.001) = 0.001;

% Flat projection
R2 = R ./ cos(PHI) .^ .2;

X_orig = R2 .* cos(TH);
Y_orig = R2 .* sin(TH);

% --- 1. Setup Interpolation Function (F) ---
% F maps (X, Y) from the original data to its W_orig value.
X = X_orig(:);
Y = Y_orig(:);
W = W_orig(:);

% Use 'natural' interpolation for smooth transitions. 'none' extrapolation
% ensures we only get values from within the data hull, returning NaN otherwise.
F = scatteredInterpolant(X, Y, W, 'natural', 'none');

% --- 2. Define Inverse Rotation Matrix ---
% To rotate the pattern +angle_deg, we must query the data at a position
% rotated by -angle_deg (inverse mapping).
angle_rad = deg2rad(-angle_deg);

R_inv = [cos(angle_rad), -sin(angle_rad);
    sin(angle_rad),  cos(angle_rad)];

% Combine original channel coordinates (Target positions)
Target_Coords = [X_orig, Y_orig]';

% --- 3. Apply Inverse Rotation and Query ---

% Source_Coords are the locations where the data was originally located
% before the pattern was conceptually rotated.
Source_Coords = R_inv * Target_Coords;

X_source = Source_Coords(1, :)';
Y_source = Source_Coords(2, :)';

% Query the original data field (F) at the new Source locations.
% W_rotated(i) gets the weight from the location that has been rotated
% into Target_Coord(i).
W_rotated = F(X_source, Y_source);

% --- 4. Handle Extrapolation (NaNs) ---
% Extrapolation is set to 'none' in Step 1, which results in NaNs for points
% outside the convex hull of the electrodes.
num_nan = sum(isnan(W_rotated));
if num_nan > 0
    % Replace NaNs with zero (a common practice, but mean or median might also be used)
    W_rotated(isnan(W_rotated)) = 0;
    fprintf('Warning: Replaced %d NaN values with 0 (extrapolated outside electrode hull).\n', num_nan);
end

fprintf('Pattern rotation complete. New weights assigned to W_rotated.\n');
end