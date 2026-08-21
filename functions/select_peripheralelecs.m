function selected_electrodes = select_peripheralelecs(EEG)
% Extract only EEG channels
eegchans = strcmp({EEG.chanlocs.type}, 'EEG');
chanLocs = EEG.chanlocs(eegchans);

x = [chanLocs(:).X]';
y = [chanLocs(:).Y]';
z = [chanLocs(:).Z]';

% Flatten 3D to 2D
z2 = z - max(z);
hypotxy = hypot(x,y);
R   = hypot(hypotxy,z2);
PHI = atan2(z2,hypotxy);
TH  = atan2(y,x);

% Remove the too small values for PHI
PHI(PHI < 0.001) = 0.001;

% Flat projection
R2 = R ./ cos(PHI) .^ .2;
X = R2.*cos(TH);
Y = R2.*sin(TH);

% Compute the geometric centre (centroid) of the 2D point cloud
centre_X = mean(X);
centre_Y = mean(Y);
% centre_X = 0;
% centre_Y = 0;

% Compute 2D distances from this empirical centre
distances_2D = sqrt((X - centre_X).^2 + (Y - centre_Y).^2);

% Find electrodes further than the radius threshold
selected_electrodes1 = distances_2D > 120;
selected_electrodes2 = Y > 100 | Y < -100; % Temporal elecs
% selected_electrodes3 = X > 60 & distances_2D > 80;

% selected_electrodes  = selected_electrodes1 | selected_electrodes2 | selected_electrodes3;
selected_electrodes  = selected_electrodes1 | selected_electrodes2;

% % Display results
% fprintf('Selected electrode indices: %s\n', mat2str(find(selected_electrodes)));
% figure; hold on;
% scatter(Y,X);
% scatter(Y(selected_electrodes),X(selected_electrodes),"filled");
% % Plot the calculated centre as a red cross for verification
% plot(centre_Y, centre_X, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
% xlim([-200 200]); ylim([-200 200]); axis square;
end