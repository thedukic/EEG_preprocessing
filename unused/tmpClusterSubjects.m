
% load('psdSpreadAll.mat');
psdSpreadAll0 = psdSpreadAll;

psdSpreadAll = psdSpreadAll';
psdSpreadAll = psdSpreadAll - min(psdSpreadAll(:));
psdSpreadAll = psdSpreadAll / max(psdSpreadAll(:));
psdSpreadAll = bsxfun(@minus, psdSpreadAll, mean(psdSpreadAll, 1));
psdSpreadAll = psdSpreadAll';
figure; histogram(psdSpreadAll(:));

% =========================================================================
[coeff, score, latent, tsquared, explained, mu] = pca(psdSpreadAll);

explained2 = explained ./ sum(explained);
Nkeep = 3;

figure; tiledlayout(1,3);
nexttile; bar(explained2); title(num2str(sum(explained2(1:Nkeep))));
nexttile; scatter(score(:,1),score(:,2));
nexttile; scatter3(score(:,1),score(:,2),score(:,3));

% =========================================================================
% NSUB = size(psdSpreadAll,1);
% nepsis = 50;
% epsis  = linspace(1,6,nepsis);
% qvec   = NaN(nepsis,1);
% 
% for epi = 1:length(epsis)
%     % scan
%     freqbands = dbscan(score(:,1:2),epsis(epi),2);
% 
%     % compute q
%     qtmp = zeros(max(freqbands),1);
%     MA = false(NSUB,1);
%     for i = 1:max(freqbands)
%         M = false(NSUB,1);
%         M(freqbands==i) = 1;
%         qtmp(i) = mean(mean(score(M,1:2))) / mean(mean(score(~M,1:2)));
%         MA = MA + M;
%     end
%     qvec(epi) = mean(qtmp) + log(mean(MA(:)));
% end
% 
% % run it again on the best epsilon value
% [~,epsiidx] = findpeaks(qvec,'NPeaks',1,'SortStr','descend');
% if isempty(epsiidx), epsiidx = round(nepsis/2); end

% =========================================================================
y = dbscan(score(:,1:Nkeep),0.25,3);

% figure;
% gscatter(score(:,1),score(:,2),y);

figure; tiledlayout(2,2);
nexttile; gscatter(score(:,1),score(:,2),y);
% nexttile; scatter(score(:,1),score(:,2));
nexttile; scatter3(score(:,1),score(:,2),score(:,3));
m = y==-1;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],num2str(sum(m)),nexttile); colorbar;
m = y==1;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],num2str(sum(m)),nexttile); colorbar;

% =========================================================================
% Step 1: Calculate the mean vector of the data
X = score(:,1:Nkeep);
mu = mean(X);

% Step 3: Calculate the Mahalanobis distance for each subject
% Using MATLAB's built-in function 'mahal'
D_M = mahal(X, X); 

% Optional: Identify potential outliers
threshold = chi2inv(0.975, size(X, 2)); % 97.5% confidence interval 
outliers = find(D_M > threshold);

disp('Potential outliers (subject indices):');
disp(subjects(outliers)');

figure; tiledlayout(2,ceil(length(outliers)/2));
for i = 1:length(outliers)
    tmp = psdSpreadAll(outliers(i),:);
    mytopoplot(tmp,[],subjects{outliers(i)},nexttile); colorbar;
end

% =========================================================================
[y,S,F,ydata,alpha] = SIMLR(psdSpreadAll,2,10,false,true);

figure; tiledlayout(2,2);
nexttile; gscatter(F(:,1),F(:,2),y);
nexttile; gscatter(ydata(:,1),ydata(:,2),y);

m = y==1;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],num2str(sum(m)),nexttile); colorbar;
m = y==2;
tmp = mean(psdSpreadAll(m,:),1);
mytopoplot(tmp,[],num2str(sum(m)),nexttile); colorbar;

% =========================================================================
% [reduction, umap, clusterIds, extras] = run_umap(psdSpreadAll);
