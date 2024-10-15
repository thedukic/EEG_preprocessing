function estimate_ictemplates(myPaths,subjects,rnum)

% IClabel = 'ECG';
IClabel = 'Blink';
% IClabel = 'Saccade';

% =========================================================================
fprintf('Estimating the %s signal reconstruction weights...\n',IClabel);
if isnumeric(rnum), rnum = num2str(rnum); end

NSUB = length(subjects);
fprintf('Loading %d datasets... Wait...\n',NSUB);

myCmap = brewermap(128,'*RdBu');
figure;

ICs = [];
cnt = 0;
for i = 1:NSUB
    load(fullfile(myPaths.preproc,subjects{i},[subjects{i} '_' myPaths.visit '_' myPaths.task '_cleandata_' rnum 'b.mat']),'EEG');

    switch IClabel
        case 'ECG'
            ICtmp = EEG.ALSUTRECHT.ica.combi.report==4;
        case 'Blink'
            ICtmp = EEG.ALSUTRECHT.ica.ICsMostLikelyBlink;
        case 'Saccade'
            ICtmp = EEG.ALSUTRECHT.ica.ICsMostLikelySaccade;
    end

    ICtmp = find(ICtmp);
    if ~isempty(ICtmp)
        theseIC = EEG.icawinv(:,ICtmp);
        ICs = [ICs, theseIC];

        for j = 1:length(ICtmp)
            cnt = cnt+1;
            nexttile;
            theseICtmp = theseIC(:,j);
            topoplot(theseICtmp,EEG.chanlocs,'maplimits',max(abs(theseICtmp))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map','shading','interp');

            title({[subjects{i} ', IC' num2str(ICtmp(j))]});
            drawnow;
        end

        fprintf('%d : %d\n',i,length(ICtmp));
    else
        fprintf('%d : %d\n',i,length(ICtmp));
    end
end

% =========================================================================
% N = vecnorm(ICs,2,1); % all the same -> 11.32
% N = vecnorm(EEG.icawinv,2,1); % all the same -> 11.32

R = corrcoef(ICs);
% figure; imagesc(R);
Rm = mean(abs(R),2);

ICbadmask = Rm<0.5;
ICs(:,ICbadmask) = [];
fprintf('Odd %s ICs removed (%d/%d)\n',IClabel,sum(ICbadmask),length(ICbadmask));

% Make sure that all the ICs are in the same direction
[s,v,d] = svd(ICs');
flipDir = sign(s(:,1));
assert(all(flipDir~=0));
ICaligned = (flipDir.*ICs')';

% Plot them again
fh = figure;
th = tiledlayout("flow");
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:size(ICaligned,2)
    nexttile;
    topoplot(ICaligned(:,i),EEG.chanlocs,'maplimits',max(abs(ICaligned(:,i)))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map','shading','interp');
    drawnow;
end
title(th,'Aligned ICs');

% The final estimate
ICweights = mean(ICaligned,2);
mytopoplot(ICweights,[],'Final reconstruction weights');

% Save
switch IClabel
    case 'ECG'
        Heartweights = ICweights;
        save(fullfile(myPaths.mycodes,'files','Heartweights'),'Heartweights');
    case 'Blink'
        Blinkweights = ICweights;
        save(fullfile(myPaths.mycodes,'files','Blinkweights'),'Blinkweights');
    case 'Saccade'
        Saccadeweights = ICweights;
        save(fullfile(myPaths.mycodes,'files','Saccadeweights'),'Saccadeweights');
end

end