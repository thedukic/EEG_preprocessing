function subject = report_proc(subject, myPaths, myTask)

fixTimeStamp = @(thisTime) strrep(strrep(char(thisTime),':','-'),' ','-');

switch myTask
    case 'open'
        % Time
        subject.time.t0 = datetime("now");

        % subject.time.procTimeTags = {myPaths.proctime; strrep(strrep(char(subject.time.t0),':','-'),' ','-')};

        % Open the report
        subject.fid = fopen(fullfile(subject.preproc, ['report_preprocess_v' myPaths.codever '.txt']), 'w+');

        % Print into the report
        fprintf(subject.fid,'%s | %s | %s dataset\n\n', myPaths.group, subject.id, myPaths.task);
        fprintf(subject.fid,'Code version %s\n', myPaths.codever);
        fprintf(subject.fid,'Started: %s\n', subject.time.t0);

    case 'close'
        % Time
        subject.time.t1 = datetime("now");
        subject.time.dd = round(minutes(diff([subject.time.t0 subject.time.t1])));

        subject.time.procTimeTags = {...
            myPaths.proctime; ...
            fixTimeStamp(subject.time.t0);
            fixTimeStamp(subject.time.t1); ...
            };

        % Print into the report
        fprintf(subject.fid,'\n\nFinished: %s\n', subject.time.t1);
        fprintf(subject.fid,'Running time: %d min.\n', subject.time.dd);

        % Close the report
        fclose(subject.fid);

        % Report
        fprintf('Finished: %s\n', subject.time.t1);
        fprintf('Running time: %d min.\n\n', subject.time.dd);

end

end