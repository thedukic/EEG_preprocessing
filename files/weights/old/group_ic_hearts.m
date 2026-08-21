
% Make sure that all the ICs are in the same direction
[s, v, d] = svd(ICs');

flipDir = sign(s(:,1));
assert(all(flipDir ~= 0));

ICaligned = (flipDir .* ICs')';

ICaligned = zscore(ICaligned);

idx = kmeans(ICaligned',3,'MaxIter',10000,'Display','final','Replicates',15);

% IChX = ICaligned(:,idx == 3);
% figure;
% for i = 1:size(IChX,2)
%     mytopoplot(IChX(:,i), [], [], nexttile); drawnow;
% end

ICh1 = mean(ICaligned(:,idx == 1),2);
ICh2 = mean(ICaligned(:,idx == 2),2);
ICh3 = mean(ICaligned(:,idx == 3),2);

figure; tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");
mytopoplot(ICh1, [], [], nexttile);
mytopoplot(ICh2, [], [], nexttile);
mytopoplot(ICh3, [], [], nexttile);