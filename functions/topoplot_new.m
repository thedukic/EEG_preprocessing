function [handle, Zi, grid_out, Xi, Yi] = topoplot_new(Values, loc_file, varargin)
% TOPOPLOT - Standalone, dependency-free mirror of the EEGLAB topoplot function.

% 1. Set Defaults
maplimits        = 'absmax';
grid_scale       = 67;
style            = 'both';
electrodes       = 'on';
headrad          = 0.5;
colormap_matrix  = parula;
plotrad_provided = false;
shading_style    = 'flat';

% Initialise
handle = []; grid_out = []; Xi = []; Yi = []; Zi = [];

% Marker Defaults
emarker_char = '.'; emarker_color = 'k'; emarker_size = 8; emarker_lw = 1;
emarker2_chans = []; emarker2_char = 'o'; emarker2_color = 'r'; emarker2_size = 10; emarker2_lw = 1;

% 2. Parse varargin
for i = 1:2:length(varargin)
    if i+1 > length(varargin); break; end
    param = lower(varargin{i});
    val = varargin{i+1};
    switch param
        case 'maplimits'
            maplimits = val;
        case 'gridscale'
            grid_scale = val;
        case 'style'
            style = lower(val);
        case 'shading'
            shading_style = lower(val);
        case 'electrodes'
            electrodes = lower(val);
        case 'plotrad'
            plotrad = val;
            plotrad_provided = true;
        case 'headrad'
            headrad = val;
        case 'colormap'
            colormap_matrix = val;
        case 'emarker'
            if iscell(val)
                if length(val)>=1 && ~isempty(val{1}); emarker_char = val{1}; end
                if length(val)>=2 && ~isempty(val{2}); emarker_color = val{2}; end
                if length(val)>=3 && ~isempty(val{3}); emarker_size = val{3}; end
                if length(val)>=4 && ~isempty(val{4}); emarker_lw = val{4}; end
            end
        case 'emarker2'
            if iscell(val)
                if length(val)>=1 && ~isempty(val{1}); emarker2_chans = val{1}; end
                if length(val)>=2 && ~isempty(val{2}); emarker2_char = val{2}; end
                if length(val)>=3 && ~isempty(val{3}); emarker2_color = val{3}; end
                if length(val)>=4 && ~isempty(val{4}); emarker2_size = val{4}; end
                if length(val)>=5 && ~isempty(val{5}); emarker2_lw = val{5}; end
            end
        case 'whitebk'
            whitebk = lower(val);
    end
end

% Check if reference channel was omitted in Values
if length(Values) == length(loc_file) - 1
    loc_file = loc_file(1:end-1);
end

% 3. Extract Coordinates from chanlocs structure
if isstruct(loc_file)
    if isfield(loc_file, 'theta') && isfield(loc_file, 'radius') && ~isempty(loc_file(1).theta)
        th = pi/180 * [loc_file.theta];
        rd = [loc_file.radius];
        x = rd .* sin(th); % X is left/right
        y = rd .* cos(th); % Y is anterior/posterior
    elseif isfield(loc_file, 'X') && isfield(loc_file, 'Y') && isfield(loc_file, 'Z')
        X = [loc_file.X]; Y = [loc_file.Y]; Z = [loc_file.Z];
        [th, el, r] = cart2sph(X, Y, Z);
        rd = pi/2 - el;
        x = rd .* cos(th);
        y = rd .* sin(th);
    else
        error('loc_file must contain theta/radius or X/Y/Z fields.');
    end
else
    error('This standalone version requires a struct (e.g., EEG.chanlocs) for loc_file.');
end

% Auto-calculate plotrad if not provided (fixes the "skirt" cutoff issue)
if ~plotrad_provided
    plotrad = max(0.5, max(rd) * 1.02);
end

% Adjust headrad if set to 'rim' so the head outline encompasses the plotrad
if ischar(headrad)
    if strcmpi(headrad, 'rim')
        headrad = plotrad;
    else
        headrad = 0; % Default fallback for unsupported string values
    end
end

Values = Values(:); x = x(:); y = y(:);

% Handle invalid data points
valid_chans = find(~isnan(Values) & ~isinf(Values));
Values = Values(valid_chans);
x = x(valid_chans);
y = y(valid_chans);

% Standardise coordinates to fit within the plotting radius
rmax = 0.5;
squeezefac = rmax / plotrad;
x = x * squeezefac;
y = y * squeezefac;

% 4. Interpolate data using MATLAB's biharmonic spline ('v4')
% Forcing symmetric limits fixes the A-P squish
xi = linspace(-rmax, rmax, grid_scale);
yi = linspace(-rmax, rmax, grid_scale);
[Xi, Yi] = meshgrid(xi, yi);

Zi = griddata(x, y, double(Values), Xi, Yi, 'v4');

% Mask data outside the plotting circle
mask_grid = (sqrt(Xi.^2 + Yi.^2) <= rmax);
Zi(~mask_grid) = NaN;
grid_out = plotrad;

% 5. Determine Colour Limits
if ischar(maplimits)
    if strcmp(maplimits, 'absmax')
        amax = max(abs(Zi(:)));
        amin = -amax;
    elseif strcmp(maplimits, 'maxmin') || strcmp(maplimits, 'minmax')
        amin = min(Zi(:));
        amax = max(Zi(:));
    else
        amin = min(Zi(:));
        amax = max(Zi(:));
    end
elseif length(maplimits) == 2
    amin = maplimits(1);
    amax = maplimits(2);
end

% 6. Render the Plot
cla; hold on;

if ~strcmp(style, 'blank')
    if strcmp(style, 'contour')
        [~, handle] = contour(Xi, Yi, Zi, 6, 'k');
    elseif strcmp(style, 'both')
        handle = surface(Xi, Yi, zeros(size(Zi))-0.1, Zi, 'EdgeColor', 'none', 'FaceColor', shading_style);
        contour(Xi, Yi, Zi, 6, 'k');
    elseif strcmp(style, 'map')
        handle = surface(Xi, Yi, zeros(size(Zi))-0.1, Zi, 'EdgeColor', 'none', 'FaceColor', shading_style);
    elseif strcmp(style, 'fill')
        [~, handle] = contourf(Xi, Yi, Zi, 6, 'k');
    end
    caxis([amin amax]);
    colormap(gca, colormap_matrix);
end

% 7. Draw Masking Ring, Head, Ears, and Nose

% --- Anti-aliasing / Masking Ring ---
% This draws a thick white ring exactly at rmax (0.5) to cleanly
% cover the staircase artifact at the edge of the data grid.
circ = linspace(0, 2*pi, 201);
rx = sin(circ); ry = cos(circ);

% Determine background color for the mask
if exist('whitebk', 'var') && strcmpi(whitebk, 'on')
    bg_color = [1 1 1];
else
    try
        bg_color = get(gcf, 'Color');
    catch
        bg_color = [1 1 1]; % Fallback to white
    end
end

% Plot the thick masking line slightly outside rmax to shave off the jaggies
% (Z-value is set to -0.05 so it sits directly on top of the surface but under the lines)
plot3(rx * (rmax + 0.015), ry * (rmax + 0.015), zeros(size(rx)) - 0.05, 'Color', bg_color, 'LineWidth', 3);
% -------------------------------------

if headrad > 0
    hwidth = 0.007;
    hin = squeezefac * headrad * (1 - hwidth/2);

    % Head outline
    plot(rx * hin, ry * hin, 'k', 'LineWidth', 2);

    % Nose
    basex = 0.18 * rmax;
    tip = 1.15 * rmax;
    tiphw = 0.04 * rmax;
    tipr = 0.01 * rmax;
    sf = headrad / plotrad;
    nose_x = [basex; tiphw; 0; -tiphw; -basex] * sf;
    nose_y = [rmax - 0.0046; tip - tipr; tip; tip - tipr; rmax - 0.0046] * sf;
    plot(nose_x, nose_y, 'k', 'LineWidth', 2);

    % Ears
    q = 0.04;
    EarX = [.492  .510  .518  .5299 .5419  .54    .547   .532   .510   .484] * sf;
    EarY = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199] * sf;
    plot(EarX, EarY, 'k', 'LineWidth', 2);
    plot(-EarX, EarY, 'k', 'LineWidth', 2);
end

% 8. Plot Electrodes (Applying emarker and emarker2 styling)
if strcmp(electrodes, 'on')
    % Intersect to find which valid channels belong to emarker2
    [~, mark2_idx] = intersect(valid_chans, emarker2_chans);
    mark1_idx = setdiff(1:length(x), mark2_idx);

    % Plot primary markers
    if ~isempty(mark1_idx)
        plot(x(mark1_idx), y(mark1_idx), emarker_char, 'Color', emarker_color, ...
            'MarkerSize', emarker_size, 'LineWidth', emarker_lw);
    end
    % Plot secondary markers (often bold/filled for highlighted subsets)
    if ~isempty(mark2_idx)
        plot(x(mark2_idx), y(mark2_idx), emarker2_char, 'Color', emarker2_color, ...
            'MarkerFaceColor', emarker2_color, 'MarkerSize', emarker2_size, 'LineWidth', emarker2_lw);
    end
elseif strcmp(electrodes, 'labels')
    if isfield(loc_file, 'labels')
        lbls = {loc_file(valid_chans).labels};
        for i = 1:length(x)
            text(x(i), y(i), lbls{i}, 'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    end
end

% 9. Finalize Axes
axis square;
axis off;
set(gca, 'Xlim', [-rmax rmax]*1.3, 'Ylim', [-rmax rmax]*1.3);
hold off;

end