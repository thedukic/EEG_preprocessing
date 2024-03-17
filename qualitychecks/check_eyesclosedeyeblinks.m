function EEG = check_eyesclosedeyeblinks(EEG)
%
% Estimate blinks in EO and compare with detections EC ???
%

% % Select these
% chaneeg  = strcmp({EEG.chanlocs.type},'EEG') & contains({EEG.chanlocs.labels},'C');
% chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
% dataeeg  = EEG.data(chaneeg,:);
% dataeog  = EEG.data(chaneog,:);
% 
% % Only EC
% ec_mask = EEG.ALSUTRECHT.blockinfo.rs_mask(~EEG.ALSUTRECHT.blockinfo.eo_mask,:);
% ec_mask = any(ec_mask,1);
% assert(length(ec_mask)==size(dataeeg,2));
% 
% % dataeeg = dataeeg(:,ec_mask);
% dataeog = dataeog(:,ec_mask);
% times   = EEG.times(ec_mask)/1000;
% times   = times-times(1);
% 
% %% Strong EMG
% % EOG: Temporarily filter, lowpass 6 Hz
% [bl, al] = butter(2,5/(EEG.srate/2),'low');
% assert(isstable(bl,al));
% dataeog = filtfilt(bl,al,dataeog);
% 
% % Detect eye blinks
% treshold = prctile(dataeog,75,2) + 3*iqr(dataeog,2);
% [qrspeaks,locs] = findpeaks(dataeog,times,'MinPeakHeight',treshold);
% 
% NEOG = length(locs);
% if NEOG>0
%     % figure; hold on;
%     % plot(times(1:100*256),dataeog(1:100*256));
%     % plot(locs(1:11),qrspeaks(1:11),'ro');
% 
%     badEpoch2 = locs*EEG.srate;
%     badEpoch2 = [badEpoch2-80; badEpoch2+80]';
%     N = length(badEpoch2(1,1):badEpoch2(1,2));
% 
%     badEpoch2(badEpoch2<1) = 1;
%     badEpoch2(badEpoch2>length(dataeog)) = length(dataeog);
% 
%     EOG = NaN(NEOG,N);
%     for i = 1:NEOG
%         if length(badEpoch2(i,1):badEpoch2(i,2))==N
%             EOG(i,:) = dataeog(badEpoch2(i,1):badEpoch2(i,2));
%         end
%     end
%     EOG(isnan(EOG(:,1)),:) = [];
%     mEOG = median(EOG);
% 
%     fh = figure; hold on;
%     F = (0:size(EOG,2)-1)./EEG.srate;
%     plot(F,EOG,'LineWidth',1.2);
%     set(gca, 'ColorOrder',brewermap(size(EOG,1),'BuGn'));
%     plot(F,mEOG,'Color',[0.8 0.1 0.1],'LineWidth',3);
%     title(['N = ' num2str(NEOG)]);
%     xlabel('Time (s)'); ylabel('VEOG amplitude');
% 
%     % Save
%     plotX=35; plotY=20;
%     set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
%     set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
%     print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_VEOG_EC']),'-dtiff','-r400');
%     close(fh);
% 
%     warning('This participant has eye blinks detected during EC!');
% end

NEOG = 0;

% Log
EEG.ALSUTRECHT.blockinfo.ec_blinks = NEOG;

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Eys-closed eye blinks\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %d\n', EEG.ALSUTRECHT.blockinfo.ec_blinks);

% data2 = dataeog;
% data2(mask)  = 0;
% data2(~mask) = 20;
%
% data2 = conv(data2,ones(1,2),'same');
% data2 = movmean(data2',16)';
%
% % EOG mask looks like it is lagging wrt EEG mask by ~d samples
% d = 20;
% data2(1:d) = [];
% data2 = [data2, data2(end)*ones(1,d)];
% assert(length(data2)==length(dataeog));
%
% % EEG: Temporarily filter, highpass 25 Hz
% [bl, al] = butter(5,5/(EEG.srate/2),'low');
% assert(isstable(bl,al));
% dataeeg = filtfilt(bl,al,dataeeg);
% dataeeg = [zeros(sum(chaneeg),1), diff(dataeeg')'];
% dataeeg = abs(dataeeg);
%
% % % Check
% % EEG0 = EEG;
% % EEG0.data(chaneeg,:) = dataeeg;
% % vis_artifacts(EEG,EEG0);
%
% treshold = prctile(dataeeg,75,2) + 2*iqr(dataeeg,2);
% mask = dataeeg<treshold;
%
% data1 = dataeeg;
% data1(mask)  = 0;
% data1(~mask) = 1;
%
% for i = 1:size(data1,1)
%     data1(i,:) = conv(data1(i,:),ones(1,10),'same');
% end
% data1 = movmean(data1',64)';
%
% % At least 25% of electrodes must be affected
% extremeMask = data1>0.1*max(data1,[],"all");
% % extremeMask = data1>0;
% extremeMask = 100*mean(extremeMask,1);
% % extremeMask(extremeMask<25) = 0;
% % extremeMask = movmean(extremeMask,64)>0;
%
% % Check
% EEG0 = EEG;
% EEG0.data = EEG0.data*0;
% EEG0.data(chaneeg,:) = 10*data1;
% EEG0.data(chaneog,:) = 10*data2;
% EEG0.data(end,:)     = 50*extremeMask;
% vis_artifacts(EEG,EEG0);

end