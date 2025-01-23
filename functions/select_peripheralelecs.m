function selected_electrodes = select_peripheralelecs(EEG)

% chanLocs = readlocs('biosemi128_eeglab.ced');
eegchans = strcmp({EEG.chanlocs.type},'EEG');
chanLocs = EEG.chanlocs(eegchans);

x = [chanLocs(:).X]';
y = [chanLocs(:).Y]';
z = [chanLocs(:).Z]';

% Flatten 3D -> 2D
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

% Compute 2D distances from the reference point
electrode_coords_2D = [Y, X];
distances_2D = sqrt(sum((electrode_coords_2D - electrode_coords_2D(1,:)).^2, 2));

% Find electrodes further than the radius threshold
selected_electrodes1 = distances_2D > 100;
selected_electrodes2 = Y>80 | Y<-80;
selected_electrodes3 = X>60 & distances_2D > 80;
selected_electrodes  = selected_electrodes1 | selected_electrodes2 | selected_electrodes3;

% % Display results
% fprintf('Selected electrode indices: %s\n', mat2str(find(selected_electrodes)));

% figure; hold on;
% scatter(Y,X); 
% scatter(Y(selected_electrodes),X(selected_electrodes),"filled"); 
% xlim([-200 200]); ylim([-200 200]);

end