function [powout, freq] = checkpowerspectrum(data,theseRows,theseTrials)

if nargin == 1
    theseRows = 1:size(data.data,1);
    theseTrials = [];
elseif nargin == 2
    theseTrials = [];
end

rsdata = false;
if isfield(data,'fsample')
    NTRL = length(data.trial);
    if NTRL==1
        rsdata = true;
    end

    fs = data.fsample;
    % data = mean(cat(3,data.trial{:}),3);
    if isempty(theseTrials)
        data = cat(3,data.trial{:});
    else
        data = cat(3,data.trial{theseTrials});
    end

else
    NTRL = size(data.data,3);
    if NTRL==1
        rsdata = true;
    end

    fs = data.srate;
    data = data.data;
    if ~isempty(theseTrials)
        data = data(:,:,theseTrials);
    end
end

NTRL = size(data,3);
if rsdata
    L = 2*fs;
    NTRL = floor(size(data,2)/(L));
    data = reshape(data(:,1:NTRL*L),size(data,1),[],NTRL);
else
    L = size(data,2);
end

if rem(L,2)>0
    L2=L+1;
else
    L2=L;
end

NFRQ = L2/2+1;
freq = fs*(0:(L2/2))/L2;

figure; hold on;
% if length(theseRows)<=5
%     tiledlayout(1,length(theseRows));
% else
%     tiledlayout(4,5);
% end

for i = 1:length(theseRows)
    pow = NaN(NTRL,NFRQ);
    for j = 1:NTRL
        tmp = fft(hanning(L)' .* squeeze(data(theseRows(i),:,j)));
        pow(j,:) = abs(tmp(1:NFRQ)).^2;
    end
    powout(i,:) = mean(pow,1);

    % nexttile;
    % plot(freq,log10(pow));
    % pow = pow./sum(pow);
    plot(freq,powout(i,:));
    title(['COMP' num2str(theseRows(i))]);
    % xlim([2 max(freq)]);
    xlim([0 100]);
    pbaspect([1.62 1 1]);

    % powall(i,:) = pow;
end

% figure;
% if length(NC)>1 && length(NC)<=5
%     tiledlayout(1,length(NC));
% else
%     tiledlayout(4,5);
% end
%
% for i = 1:length(NC)
%     nexttile;
%     plot(data(NC(i),:));
%     title(['COMP' num2str(NC(i))]);
%     pbaspect([1.62 1 1]); axis tight;
% end

end