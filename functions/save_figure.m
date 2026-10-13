function save_figure(fh, plot_dir, plot_name, plot_dims, plot_res)
% SAVE_FIGURE Exports figure at specified DPI and closes it.
%
% Inputs:
%   fh        : Figure handle
%   save_dir  : Output directory path
%   file_name : Output filename (without extension)
%   plot_dims : [width, height] in centimetres (optional)
%   res_dpi   : Output resolution in DPI (optional, default: 200)

if nargin < 5 || isempty(plot_res)
    plot_res = 200;
end

if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

% Set canvas properties
fh.Color = 'w';

% Adjust physical dimensions if explicitly passed
if nargin >= 4 && ~isempty(plot_dims)
    fh.Units = 'centimeters';
    fh.Position(3) = plot_dims(1); % Width  (43.5 cm)
    fh.Position(4) = plot_dims(2); % Height (15.0 cm)
end

% Ensure layout engines finish rendering
drawnow;

% Export to TIFF
full_path = fullfile(plot_dir, [plot_name '.tif']);
exportgraphics(fh, full_path, 'Resolution', plot_res);

% Delete figure handle to release memory
delete(fh);

end

% function save_figure(fh, save_path, file_name, plot_dims, resolution)
% % SAVE_AND_CLEAR_FIGURE Prints an un-referenced or hidden figure to disk
% % and aggressively forces MATLAB and Java to release graphics memory.
% %
% % Inputs:
% %   fh        : Figure handle (e.g., fh)
% %   save_path : Absolute path to target folder (string)
% %   file_name : Output filename without extension (string)
% %   plot_dims : 1x2 numeric vector [width, height] in centimeters
% %   resolution: Numeric scalar for DPI (default = 150)
%
% if nargin < 5 || isempty(resolution)
%     resolution = 150;
% end
%
% % 1. Enforce precise output geometry
% plotX = plot_dims(1);
% plotY = plot_dims(2);
%
% set(fh, 'InvertHardCopy', 'Off', 'Color', [1 1 1]);
% set(fh, 'PaperPositionMode', 'Manual', ...
%     'PaperUnits', 'Centimeters', ...
%     'PaperPosition', [0 0 plotX plotY], ...
%     'PaperSize', [plotX plotY]);
%
% % 2. Ensure target path exists
% if ~exist(save_path, 'dir')
%     mkdir(save_path);
% end
%
% % 3. Print file to disk using high-quality TIFF
% full_output_file = fullfile(save_path, file_name);
% print(fh, full_output_file, '-dtiff', sprintf('-r%d', resolution));
%
% % 4. Hard purge from memory
% % Bypasses the UI close queue to instantly destroy the graphics object
% delete(fh);
%
% % Flush graphics queue and force Java Garbage Collection
% drawnow;
% java.lang.System.gc();
% end