function EEG = do_jointdecorrelation(EEG,EXT,cfg)

NREMOVE = 1;

% sr = EEG.srate;
maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
chanveog = find(strcmp({EXT.chanlocs.labels},'VEOG'));
chanheog = find(strcmp({EXT.chanlocs.labels},'HEOG'));
VEOG = EXT.data(chanveog,:);
% HEOG = EXT.data(chanheog,:);

% [B,A] = butter(2,1/(sr/2), 'high');
% tmp = filter(B,A,VEOG);
%
% treshold = 8 * median(abs(tmp));
% mask = abs(tmp) > treshold;

% eye_channels=[93 94 81 82 71 72];
% VEOG = x(:,eye_channels);
% tmp = nt_pca(filter(B,A,VEOG));
% mask = abs(tmp(:,1))>3*median(abs(tmp(:,1)));

% Temporarily filter for better detection
[bl, al] = butter(2,5/(EEG.srate/2),'low');
% [bh, ah] = butter(2,1/(EEG.srate/2),'high');

eyeBlinkVEOG = do_filteringcore(bl,al,VEOG,EEG.event,EEG.srate);
% eyeBlinkVEOG = do_filteringcore(bh,ah,eyeBlinkVEOG,EEG.event,EEG.srate);
treshold = 5 * median(abs(eyeBlinkVEOG));
maskBlink1 = abs(eyeBlinkVEOG) > treshold;

if EEG.ALSUTRECHT.leftovers.flagCorrTstat
    mask1 = EEG.ALSUTRECHT.leftovers.blinksTstat > 3;
else
    mask1 = false(1,EEG.nbchan);
end
if EEG.ALSUTRECHT.leftovers.flagCorrPercent
    mask2 = EEG.ALSUTRECHT.leftovers.blinksRatio > 0.1;
else
    mask2 = false(1,EEG.nbchan);
end
maskBlink2 = mask1 | mask2;

% eyeBlinkEEG0 = EEG.data';
% eyeBlinkEEG1 = EEG.data';

% eyeBlinkEEG0 = EEG.data';
eyeBlinkEEG = do_filteringcore(bl,al,EEG.data,EEG.event,EEG.srate);
% eyeBlinkEEG = do_filteringcore(bh,ah,eyeBlinkEEG,EEG.event,EEG.srate);
eyeBlinkEEG = mean(eyeBlinkEEG(maskBlink2,:),1)';

treshold = 5 * median(abs(eyeBlinkEEG));
maskBlink2 = abs(eyeBlinkEEG) > treshold;

maskBlink = maskBlink1(:) & maskBlink2(:);

% eyeBlinkEEG1 = do_filteringcore(bl,al,EEG.data,EEG.event,EEG.srate);
% % eyeBlinkEEG1 = do_filteringcore(bh,ah,eyeBlinkEEG1,EEG.event,EEG.srate);
% eyeBlinkEEG1 = eyeBlinkEEG1';
% eyeBlinkEEG0 = eyeBlinkEEG1;

figure; hold on;
plot(eyeBlinkVEOG);
% plot(mean(eyeBlinkEEG(:,[80 81 93]),2));
plot(eyeBlinkEEG);
scatter(1:length(maskBlink),treshold*maskBlink);
axis tight;

C0 = nt_cov(eyeBlinkEEG0);
C1 = nt_cov(bsxfun(@times, eyeBlinkEEG1, double(maskBlink)'));
[todss, pwr0, pwr1] = nt_dss0(C0,C1,[],0.05);
p1 = pwr1 ./ pwr0;

figure; plot(p1, '.-');
ylabel('score'); xlabel('component'); title ('eyeblink DSS');

figure;
th = tiledlayout(2,5);
th.TileSpacing = 'compact'; th.Padding = 'compact';
for i_cmp = 1:min(10,length(p1))
    mytopoplot(todss(:,i_cmp),[],num2str(round(p1(i_cmp),3)),nexttile);
end

eye_components = eyeBlinkEEG0 * todss;
eyeBlinkEEG00 = nt_tsr(eyeBlinkEEG0, eye_components(:,1:NREMOVE));

time = 1:50*256;
figure; clf;
subplot 311; plot(eye_components(time,1:NREMOVE)); title('eye components'); axis tight;
subplot 312; plot(mean(eyeBlinkEEG0(time,maskChanBlink),2)); title('eyeblinks before removing'); axis tight;
subplot 313; plot(mean(eyeBlinkEEG00(time,maskChanBlink),2)); title('eyeblinks removed'); axis tight;

% Return
EEG.data = eyeBlinkEEG00';

end