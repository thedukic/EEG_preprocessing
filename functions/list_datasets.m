function [blocklist, NBLK] = list_datasets(datapath,thistask)

if strcmpi(thistask,'RS')
    % Resting-state
    tmp1 = dir([datapath '\*EO*']);
    tmp2 = dir([datapath '\*EC*']);
    dataname  = {tmp1.name, tmp2.name};
elseif strcmpi(thistask,'MT')
    % Motor tasks
    tmp1 = dir([datapath '\*MT2*']);
    tmp2 = dir([datapath '\*MT3*']);
    tmp3 = dir([datapath '\*MT5*']);
    dataname  = {tmp1.name, tmp2.name, tmp3.name};
else
    % SART/MMN/EO/EC
    dataname = dir([datapath '\*' thistask '*']);
    dataname = {dataname.name};
end

% Maybe not needed anymore
dataname = dataname(cellfun(@isempty,regexp(dataname,'bad')));

blocklist = fullfile(datapath,dataname);
NBLK = length(blocklist);

end