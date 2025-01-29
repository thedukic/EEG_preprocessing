
DSR   = 100; % determines granularity (& minimum cluster size)
FLAGS = []; % 'norm' or 'norm2'

x = EEG.data(1:128,:)';
nt_cluster_jd(x,DSR);
[IDX,TODSS,SCORE] = nt_cluster_jd(x,DSR,FLAGS);
% disp(['score: ',num2str(SCORE')]);

nt_bias_cluster(x,DSR)

[M,y,score] = nt_sca(x,50);
th = tiledlayout(2,5);
for i = 1:10
    mytopoplot(M(:,i),[],'',nexttile); % colorbar;
end

% =========================================================================
% eog = squeeze(EEG.data(130,:,:));
% figure; plot(EEG.times,eog);

xx = permute(EEG.data(1:128,:,:),[2 1 3]);

% DSS to emphasize repeatablity
[todss, pwr0, pwr1] = nt_dss1(xx);
fromdss = pinv(todss);
z = nt_mmat(xx,todss);

th = tiledlayout(2,5);
for i = 1:10
    mytopoplot(fromdss(i,:),[],'',nexttile); % colorbar;
end

figure; clf;
subplot 121;
plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title ('repeatability DSS');
subplot 122;
nt_bsplot(z(:,1,:)); title('best DSS component');

% denoise by selecting best components and projecting back to sensor space
NKEEP = 7;
y = nt_mmat(xx,todss(:,1:NKEEP)*fromdss(1:NKEEP,:));

figure; clf;
subplot 121; plot(mean(xx,3)); title('before DSS denoising'); axis tight; ylim([-20 15]);
subplot 122; plot(mean(y,3)); title('after'); axis tight; ylim([-20 15]);

% =========================================================================
[Y,excentricity,removed,cov1,cov2] = nt_lsp(xx,3);

figure;
subplot 131; plot(mean(xx,3)); title('before DSS denoising'); ylim([-20 15]);
subplot 132; plot(mean(Y,3)); title('after'); ylim([-20 15]);
subplot 133; plot(mean(xx,3)-mean(Y,3)); title('after'); ylim([-5 5]);

