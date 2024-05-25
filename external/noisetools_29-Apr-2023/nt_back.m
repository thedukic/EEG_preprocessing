function [iBest,y,topo]=nt_back(x,z,layout)
%[iBest,y,topo]=nt_back(x,z,layout) - back-project component
%
%  iBest: best channel per component 
%  y: data from best channel per component 
%  topoB: amplitude of other channels relative to best 
%
%  x: raw data (timeXchannels pr timeXchannelsXtrials)
%  z: component (timeX1 or timeX1Xtrials)
%  layout: for optional topoplot
% 
% NoiseTools
%

if nargin==0; test_code; return; end

if nargin<2; error('!'); end
if nargin<3||isempty(layout); layout=[]; end
if ndims(x)~=ndims(z); error('!'); end

if ndims(x)==3
    % 3D: unfold to 2D
    [m,n,o]=size(x);
    x=nt_unfold(x);
    z=nt_unfold(z);
    if nargout~=0
        [iBest,y,topo]=nt_back(x,z,layout);
        y=nt_fold(y,m);
    else
        nt_back(x,z,layout);
    end
    return;
end


% correlation with each channel
z=nt_normcol(z);
topo=nt_xcov(z,x)/size(x,1); % non-normalized

% best channel for each component
[~,iBest]=max(abs(topo)');

% back-project each component to its best channel
for iComp=1:size(z,2)
    y(:,iComp)=z(:,iComp)*topo(iComp,iBest(iComp));
end

% make sure topo has right sign 
for iComp=1:size(z,2)
    topo(iComp,:)=topo(iComp,:)*sign(mean(z(:,iComp).*y(:,iComp)));
end




if nargout==0
    disp('nt_back: no output requested, plot back matrices');
    
    J=floor(sqrt(size(z,2)));
    K=ceil(size(z,2)/J);

    figure(101); clf;
    nt_banner('best channel for each component','fontsize',16)
    for iComp=1:size(z,2)
        subplot(J,K,iComp);
        plot(x(:,iBest(iComp)), ':k'); hold on
        plot(y(:,iComp), 'b');
        title(['c',num2str(iComp), ', ch', num2str(iBest)]);
    end
    
    if ~isempty(layout)
        figure(102); clf
        nt_banner('RMS per channel','fontsize',16)
        for iComp=1:size(z,2)
            subplot(J,K,iComp);
            nt_topoplot(layout, topo(iComp,:)); colorbar
            title(iComp);
        end
    end
    
    clear iBest y topo
 end



end % function




% test/example code
function test_code
    disp('nt_back test code');
    disp('3D random matrix, apply nt_dss1, give to nt_back()');
    
    x=randn(1000,10,10);
    [todss]=nt_dss1(x);
    nt_back(x,todss)
end % function