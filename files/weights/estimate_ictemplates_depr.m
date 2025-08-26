function estimate_ictemplates(myPaths,subjects)

% Estimate IC templates for these artifacts
IClabel = {'ECG','Blink','Saccade'};

% =========================================================================
myCmap = brewermap(128,'*RdBu');

NCMP = length(IClabel);
NSUB = length(subjects);
fprintf('Loading %d datasets. Wait...\n',NSUB);

for i_type = 1:NCMP
    fprintf('Estimating the %s signal reconstruction weights...\n',IClabel{i_type});
    ICs = [];
    cnt = 0;

    % Plot
    fh = figure;
    th = tiledlayout("flow");
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    for i_sub = 1:NSUB
        % Loads
        load(fullfile(myPaths.preproc,subjects{i_sub},[subjects{i_sub} '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'b.mat']),'EEG');

        switch i_type
            case 1
                % 'ECG'
                ICtmp = EEG.ALSUTRECHT.ica.combi.report == 4;
            case 2
                % 'Blink'
                ICtmp = EEG.ALSUTRECHT.ica.ICsMostLikelyBlink;
            case 3
                % 'Saccade'
                ICtmp = EEG.ALSUTRECHT.ica.ICsMostLikelySaccade;
        end

        ICtmp = find(ICtmp);
        if ~isempty(ICtmp)
            theseIC = EEG.icawinv(:,ICtmp);
            ICs = [ICs, theseIC];

            for i_ics = 1:length(ICtmp)
                cnt = cnt+1;
                nexttile;
                theseICtmp = theseIC(:,i_ics);
                topoplot(theseICtmp,EEG.chanlocs,'maplimits',max(abs(theseICtmp))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map','shading','interp');

                title({[subjects{i_sub} ', IC' num2str(ICtmp(i_ics))]});
                drawnow;
            end

            fprintf('%d : %d\n',i_sub,length(ICtmp));
        else
            fprintf('%d : %d\n',i_sub,length(ICtmp));
        end
    end

    % =========================================================================
    % N = vecnorm(ICs,2,1); % all the same -> 11.32
    % N = vecnorm(EEG.icawinv,2,1); % all the same -> 11.32

    R = corrcoef(ICs);
    % figure; imagesc(R);
    Rm = mean(abs(R),2);

    ICbadmask = Rm < 0.5;
    ICs(:,ICbadmask) = [];
    fprintf('Odd %s ICs removed (%d/%d)\n',IClabel{i_type},sum(ICbadmask),length(ICbadmask));

    % Make sure that all the ICs are in the same direction
    [s,v,d] = svd(ICs');
    flipDir = sign(s(:,1));
    assert(all(flipDir ~= 0));
    ICaligned = (flipDir.*ICs')';

    % Plot them again
    fh = figure;
    th = tiledlayout("flow");
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    for i_sub = 1:size(ICaligned,2)
        nexttile;
        topoplot(ICaligned(:,i_sub),EEG.chanlocs,'maplimits',max(abs(ICaligned(:,i_sub)))*[-1 1],'headrad','rim','colormap',myCmap,'whitebk','on','style','map','shading','interp');
        drawnow;
    end
    title(th,[IClabel{i_type} ': Aligned ICs']);

    % The final estimate
    ICweights = mean(ICaligned,2);
    mytopoplot(ICweights,[],[IClabel{i_type} ': Final reconstruction weights']);

    % % Save
    % switch i_type
    %     case 1
    %         % 'ECG'
    %         Heartweights = ICweights;
    %         save(fullfile(myPaths.mycodes,'files','Heartweights'),'Heartweights');
    %     case 2
    %         % 'Blink'
    %         Blinkweights = ICweights;
    %         save(fullfile(myPaths.mycodes,'files','Blinkweights'),'Blinkweights');
    %     case 3
    %         % 'Saccade'
    %         Saccadeweights = ICweights;
    %         save(fullfile(myPaths.mycodes,'files','Saccadeweights'),'Saccadeweights');
    % end
end

end