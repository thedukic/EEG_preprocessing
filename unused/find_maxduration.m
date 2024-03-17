function maxdur = find_maxduration(mask)
% mask = abs(data)>treshold;

[N, T] = size(mask);
if T==1
    mask = mask';
    [N, T] = size(mask);
end

mask = [false(N,1), mask, false(N,1)];

maxdur = zeros(N,1);
for i = 1:N
    jump = find(diff(mask(i,:))~=0);
    if ~isempty(jump)
        % jump = [0, tmp, NMSK];
        % maxdur(i) = max(jump(1:2:end)-jump(2:2:end));
        maxdur(i) = max(jump(2:2:end)-jump(1:2:end));
    end
end
% maxdur = max(maxdur);

end