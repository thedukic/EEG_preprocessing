function [data, report] = interpolate_epochs(data, elec, bad_elecxtrial, ignore_chans, max_badelec)

% Assert the right format
num_trials = size(data, 3);
assert(length(bad_elecxtrial) == num_trials);

% Handle optional ignore_chans argument
if nargin < 4 || isempty(ignore_chans)
    ignore_chans = [];
end

% EEGLAB vs FieldTrip elecs definition
if isfield(elec, 'elecpos')
    assert(size(data, 1) == length(elec.label));
    isFieldTrip = true;
else
    assert(size(data, 1) == length(elec));
    isFieldTrip = false;
end

% Reformat channel positions
if isFieldTrip
    xelec = elec.elecpos(:,1)';
    yelec = elec.elecpos(:,2)';
    zelec = elec.elecpos(:,3)';
else
    xelec = [elec(:).X];
    yelec = [elec(:).Y];
    zelec = [elec(:).Z];
end

% Project to unit sphere
rad = sqrt(xelec.^2 + yelec.^2 + zelec.^2);
xelec = xelec ./ rad;
yelec = yelec ./ rad;
zelec = zelec ./ rad;

% Precompute the full G matrix for all channels once
G_full = computeg(xelec, yelec, zelec, xelec, yelec, zelec);

% Allocate
listFixed  = false(num_trials, 1);
listRemove = false(num_trials, 1);

for i = 1:num_trials
    if ~isempty(bad_elecxtrial{i})
        if length(bad_elecxtrial{i}) <= max_badelec
            % The trial can be fixed
            badchans  = bad_elecxtrial{i};

            % Exclude both bad channels and channels to ignore from the good list
            goodchans = setdiff(1:size(data, 1), [badchans, ignore_chans]);

            % Extract the relevant precomputed G matrices
            Gelec = G_full(goodchans, goodchans);
            Gsph  = G_full(badchans, goodchans);

            data(badchans, :, i) = spheric_spline(Gelec, Gsph, data(goodchans, :, i));

            listFixed(i) = true;
        else
            % The trial is too noisy
            listRemove(i) = true;
        end
    end
end

% Check
assert(~any(listRemove + listFixed == 2));

% Report
report = [];
report.listFixed  = find(listFixed);
report.listRemove = find(listRemove);

fprintf('Number of trials fixed: %d\n', sum(listFixed));
fprintf('Number of trials to be removed: %d\n', sum(listRemove));

end

% =========================================================================
% Helper functions
% =========================================================================

function allres = spheric_spline(Gelec, Gsph, values)
numpoints = size(values, 2);

% Compute solution for parameters C
meanvalues = mean(values, 1);

% MATLAB implicit expansion handles the subtraction directly
values = values - meanvalues;

values = [values; zeros(1, numpoints)];
C = pinv([Gelec; ones(1, size(Gelec, 2))]) * values;

% Vectorised application of results
allres = Gsph * C;

% Implicit expansion for adding the mean back
allres = allres + meanvalues;
end

function g = computeg(x, y, z, xelec, yelec, zelec)
% Ensure all input vectors are strictly normalised to a unit sphere (radius = 1)
r_xyz = sqrt(x.^2 + y.^2 + z.^2);
x = x ./ r_xyz; y = y ./ r_xyz; z = z ./ r_xyz;

r_elec = sqrt(xelec.^2 + yelec.^2 + zelec.^2);
xelec = xelec ./ r_elec; yelec = yelec ./ r_elec; zelec = zelec ./ r_elec;

% Compute the cosine of the angle between points
EI = x(:)*xelec(:)' + y(:)*yelec(:)' + z(:)*zelec(:)';

% Force-clamp the matrix to eliminate floating-point rounding errors
EI = max(min(real(EI), 1), -1);

% Preallocate using the explicit matrix dimensions
g = zeros(size(EI));

m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n, EI);

    % Reshape ensures safety against matrix dimension collapse
    g = g + ((2*n+1)/(n^m*(n+1)^m)) * reshape(L(1,:,:), size(EI));
end
g = g / (4*pi);
end