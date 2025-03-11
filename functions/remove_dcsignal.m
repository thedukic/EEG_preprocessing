function data = remove_dcsignal(data,windowsam,chanArray)
% ERPLAB toolbox function
% Removes mean of data (DC offset)
% Input data dimensions have to be channels x samples

if nargin < 3
    chanArray = 1:size(data,1);
end
if nargin < 2
    windowsam = [1 size(data,2)];
end
if length(windowsam) ~= 2
    windowsam = [1 windowsam(1)];
end

% Control point
if diff(windowsam) > size(data,2)
    windowsam = [1 size(data,2)];
end

% Offset
meanValue = mean(data(chanArray,windowsam(1):windowsam(2)), 2);
data = data - meanValue;

end