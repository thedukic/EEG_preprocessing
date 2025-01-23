function [data, report] = interpolate_epochs(data,elec,badElecsPerTrial,ignore_chans,maxBadChans)

% Assert the right format
assert(size(data,3) == length(badElecsPerTrial));

% EEGLAB vs FieldTrip elecs definition
if isfield(elec,'elecpos')
    assert(size(data,1) == length(elec.label));
    isFieldTrip = true;
else
    assert(size(data,1) == length(elec));
    isFieldTrip = false;
end

listFixed = [];
listNotFixed = [];

for i = 1:length(badElecsPerTrial)
    if ~isempty(badElecsPerTrial{i})
        % The trial is noisy
        if length(badElecsPerTrial{i}) <= maxBadChans
            % The trial can be fixed
            badchans  = badElecsPerTrial{i};
            goodchans = setdiff(1:size(data,1), badchans);
            % goodchans = setdiff(goodchans, ignore_chans);

            if isFieldTrip
                xelec = elec.elecpos(goodchans,1)';
                yelec = elec.elecpos(goodchans,2)';
                zelec = elec.elecpos(goodchans,3)';
                xbad  = elec.elecpos(badchans,1)';
                ybad  = elec.elecpos(badchans,2)';
                zbad  = elec.elecpos(badchans,3)';
            else
                xelec = [elec(goodchans).X];
                yelec = [elec(goodchans).Y];
                zelec = [elec(goodchans).Z];
                xbad  = [elec(badchans).X];
                ybad  = [elec(badchans).Y];
                zbad  = [elec(badchans).Z];
            end

            rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
            xelec = xelec./rad;
            yelec = yelec./rad;
            zelec = zelec./rad;

            rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
            xbad = xbad./rad;
            ybad = ybad./rad;
            zbad = zbad./rad;

            data(badchans,:,i) = spheric_spline(xelec, yelec, zelec, xbad, ybad, zbad, data(goodchans,:,i));
            listFixed = [listFixed  i];
            % fprintf('Trial %d fixed.\n',i);
        else
            % The trial is too noisy
            % fprintf('Trial %d is too noisy and thus not fixed.\n',i);
            listNotFixed = [listNotFixed i];
        end
    end
end

% Report
report = [];
report.listFixed = listFixed;
report.listNotFixed = listNotFixed;

fprintf('Number of trials fixed: %d\n',length(listFixed));
fprintf('Number of trials to be removed: %d\n',length(listNotFixed));

% =========================================================================
% =========================================================================
% Helper functions
% =========================================================================
% =========================================================================
function allres = spheric_spline(xelec, yelec, zelec, xbad, ybad, zbad, values)

newchans = length(xbad);
numpoints = size(values,2);

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);

% compute solution for parameters C
% ---------------------------------
meanvalues = mean(values);
values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero

values = [values;zeros(1,numpoints)];
C = pinv([Gelec;ones(1,length(Gelec))]) * values;
clear values;
allres = zeros(newchans, numpoints);

% apply results
for j = 1:size(Gsph,1)
    allres(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));
end
allres = allres + repmat(meanvalues, [size(allres,1) 1]);

% =========================================================================
% compute G function
% =========================================================================
function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +...
    (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
    (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));

m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n,EI);
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);
