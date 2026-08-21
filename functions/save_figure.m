function save_figure(fh, save_dir, file_name, plot_dims)
% SAVE_FIGURE Exports a figure using the modern graphics pipeline and closes it.
%
% Inputs:
%   fh        : Figure handle
%   save_dir  : Directory path to save the file
%   file_name : Name of the file (without extension)
%   plot_dims : [width, height] in centimetres (optional)

% Ensure background is white
fh.Color = 'w';

% Apply target physical dimensions before exporting
if nargin >= 4 && ~isempty(plot_dims)
    fh.Units = 'centimeters';
    fh.Position(3:4) = plot_dims;
end

% Force MATLAB to render the tiled layout geometries (now using CPU only)
drawnow;

% Construct path and export safely
full_path = fullfile(save_dir, [file_name '.tif']);
exportgraphics(fh, full_path, 'Resolution', 200);

% Destroy the figure to prevent memory leaks
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