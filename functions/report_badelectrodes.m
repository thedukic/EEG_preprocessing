function INPUTDATA = report_badelectrodes(INPUTDATA,thisType)

if strcmpi(thisType,'individual')
    assert(isstruct(INPUTDATA));

    % Log
    INPUTDATA.ALSUTRECHT.badchaninfo.badElectrodes = sort([INPUTDATA.ALSUTRECHT.badchaninfo.flatElectrodes, INPUTDATA.ALSUTRECHT.badchaninfo.PREPElectrodes]);
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'Bad  electrodes\n');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
    str = strjoin(INPUTDATA.ALSUTRECHT.badchaninfo.flatElectrodes,', ');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'Flat electrodes: %s\n', str);
    str = strjoin(INPUTDATA.ALSUTRECHT.badchaninfo.wica.fixed,', ');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'wICA electrodes: %s\n', str);
    str = arrayfun(@(x) num2str(x,'%1.1f'),INPUTDATA.ALSUTRECHT.badchaninfo.wica.pvec(INPUTDATA.ALSUTRECHT.badchaninfo.wica.ics),'uni',0);
    str = strjoin(str,', ');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'wICA P-values:   %s\n', str);
    str = strjoin(INPUTDATA.ALSUTRECHT.badchaninfo.PREPElectrodes,', ');
    fprintf(INPUTDATA.ALSUTRECHT.subject.fid,'PREP electrodes: %s\n', str);

    % Plot
    fh = figure;
    th = tiledlayout(1,3);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    % EEG electrode colours/labels
    myCmap = brewermap(128,'RdPu');
    chanlabseeg = {INPUTDATA.allchans(strcmp({INPUTDATA.allchans.type},'EEG')).labels};
    chanlocseeg = INPUTDATA.allchans(strcmp({INPUTDATA.allchans.type},'EEG'));

    mask = double(ismember(chanlabseeg,INPUTDATA.ALSUTRECHT.badchaninfo.flatElectrodes));
    nexttile;
    topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','off','style','map');
    title(['Flat, N = ' num2str(sum(mask))]); axis tight;

    mask = double(ismember(chanlabseeg,INPUTDATA.ALSUTRECHT.badchaninfo.wica.fixed));
    nexttile;
    topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','off','style','map');
    title(['wICA, N = ' num2str(sum(mask))]); axis tight;

    mask = double(ismember(chanlabseeg,INPUTDATA.ALSUTRECHT.badchaninfo.PREPElectrodes));
    nexttile;
    topoplot(mask,chanlocseeg,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','off','style','map');
    title(['PREP, N = ' num2str(sum(mask))]); axis tight;

    % Save
    plotX=20; plotY=8;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(INPUTDATA.ALSUTRECHT.subject.preproc,[INPUTDATA.ALSUTRECHT.subject.id '_badelectrodes']),'-dtiff','-r400');
    close(fh);

elseif strcmpi(thisType,'group')

    elec = biosemi128(58,true);
    maskelec = zeros(length(elec.label),1);
    rnum = 1; % HARDCODED FOR NOW 

    subjects = list_subjects(myfolders.preproc,[]);
    NSUB = length(subjects);
    N = NaN(NSUB,1);
    for i = 1:NSUB
        load([fullfile(INPUTDATA.preproc, subjects{i}) [subjects{i} '_' INPUTDATA.visit '_' INPUTDATA.task '_cleandata_' rnum '.mat']],'EEG');

        maskelec = maskelec + double(ismember(elec.label,EEG.ALSUTRECHT.badchaninfo.badElectrodes));
        N(i)= length(INPUTDATA.badElectrodes);
    end
    maskelec = maskelec./NSUB;

    % Plot
    fh = figure;
    myCmap = brewermap(128,'RdPu');
    chanlocs = readlocs('biosemi128_eeglab.ced');
    topoplot(maskelec,chanlocs,'maplimits',[0 1],'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','off','style','map');
    % title(['Flat, N = ' num2str(sum(mask))]); axis tight;
end

end