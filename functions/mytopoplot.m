function mytopoplot(data,mask,myTitle,ah,myClim)

% chanlocs = readlocs('biosemi128_eeglab.ced');
load('biosemi128_eeglab.mat','chanlocs');
% fcap = 0.5;
fcap = 'rim';

if iscell(myTitle)
    myTitlePlot      = myTitle{1};
    myTitleColourbar = myTitle{2};
else
    myTitlePlot      = myTitle;
    myTitleColourbar = '';
end

% Determine where to plot
if isempty(ah)
    figure; ah = gca;
else
    axes(ah); % gca;
end

% Colour limits
if ~exist('myClim','var')
    myDlim = [min(data), max(data)];
    myClim = max(abs(data))*[-1 1];

    % y = prctile(data,70);
    % mask = data>y;

    if diff(myDlim) == 0
        myClim(1) = 0;
    elseif all(myDlim>=0)
        % myClim(1) = 0;
        myClim(1) = myDlim(1);
    elseif all(myDlim<=0)
        % myClim(2) = 0;
        myClim(2) = myDlim(2);
    end
end
myClim = 0.98*myClim;

% Colour map
if all(data>=0)
    % myCmap = brewermap([],'RdPU');
    myCmap = brewermap([],'Reds');
elseif all(data<=0)
    myCmap = brewermap([],'*Blues');
else
    myCmap = brewermap([],'*RdBu');
end

% Example
% topoplot(dataPlot, chanlocs, 'style', 'map', 'gridscale', 300, 'maplimits', 0.95*max(abs(dataPlot))*[-1 1], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{dataMask,'o','k',4,1},'colormap',myCmap);

% Plot
if isempty(mask)
    topoplot(data,chanlocs,'headrad',fcap,'whitebk','on','electrodes','on','style','both','shading','interp','gridscale',300);
else
    if length(mask)==128
        mask = find(mask);
    end
    topoplot(data,chanlocs,'headrad',fcap,'whitebk','on','electrodes','on','style','both','shading','interp','gridscale',300,'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{mask,'o','k',4,1}); % 'emarker2',{mask,'d','k',10,1}
end

% Title
if ~isempty(myTitlePlot) || strcmpi(myTitlePlot,'')
    title(ah,myTitlePlot);
end

axis tight;
colormap(ah,myCmap);
clim(ah,myClim);

% % Colourbar
% cbh = colorbar;
% if ~isempty(myTitleColourbar)
%     cbh.Label.String = myTitleColourbar;
% end

end