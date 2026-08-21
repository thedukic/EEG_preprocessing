function mytopoplot(data, mask, myTitle, ah, myClim)

% chanlocs = readlocs('biosemi128_eeglab.ced');
load('biosemi128_eeglab.mat', 'chanlocs');

% Topoplot style
% fcap = 0.5;
fcap = 'rim';

if ~exist('myTitle','var')
    myTitle = '';
end
if ~exist('mask','var')
    mask = [];
end

% if iscell(myTitle)
%     myTitlePlot      = myTitle{1};
%     myTitleColourbar = myTitle{2};
% else
%     myTitlePlot      = myTitle;
%     myTitleColourbar = '';
% end
myTitlePlot = myTitle;

% Determine where to plot without disrupting window visibility
if exist('ah','var') && ~isempty(ah)
    % Find the parent figure of the axes and set it as active silently
    set(ancestor(ah, 'figure'), 'CurrentAxes', ah);
else
    % Create a new figure only if no valid axes handle was provided
    figure;
    ah = gca;
end

% Colour limits
if ~exist('myClim','var')
    myDlim = [min(data), max(data)];
    myClim = max(abs(data)) * [-1 1];

    % y = prctile(data,70);
    % mask = data>y;

    if diff(myDlim) == 0
        if all(data == 0)
            myClim = [0 1];
        else
            myClim(1) = 0;
        end
    elseif all(myDlim>=0)
        % myClim(1) = 0;
        myClim(1) = myDlim(1);
    elseif all(myDlim<=0)
        % myClim(2) = 0;
        myClim(2) = myDlim(2);
    end
end

myClim = 0.98 * myClim;

% Colour map
if all(data >= 0)
    % myCmap = brewermap([],'RdPU');
    myCmap = brewermap([], 'Reds');
elseif all(data<=0)
    myCmap = brewermap([], '*Blues');
else
    myCmap = brewermap([], '*RdBu');
end

% Plot
if isempty(mask)
    topoplot_new(data, chanlocs, ...
        'headrad', fcap, 'whitebk', 'on', 'electrodes', 'off', 'style', 'both', 'shading', 'interp', 'gridscale', 300, 'maplimits', myClim);
else
    if length(mask) == 128
        mask = find(mask);
    end
    topoplot_new(data, chanlocs, ...
        'headrad', fcap, 'whitebk', 'on', 'electrodes', 'on', 'style', 'both', 'shading', 'interp', 'gridscale', 300, 'maplimits', myClim, ...
        'emarker', {'.',[.5 .5 .5],[],2}, 'emarker2', {mask,'o','k',4,1});
    % 'emarker2',{mask,'d','k',10,1}
end

% Title
if ~isempty(myTitlePlot) || strcmpi(myTitlePlot,'')
    title(ah,myTitlePlot);
end

axis tight;
colormap(ah, myCmap);
clim(ah, myClim);

% % Colourbar
% cbh = colorbar;
% if ~isempty(myTitleColourbar)
%     cbh.Label.String = myTitleColourbar;
% end

end