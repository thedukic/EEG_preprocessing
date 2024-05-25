function h=nt_banner(text,varargin)
%h=nt_banner(text,varargin) - annotate with text at head of figure
% 
%  h: handle to annotation
%  
%  text: string to print
%  varargin: arguments to pass to annotation()
% 
%  default horizontal position is 'center'
% 
% NoiseTools

set(gcf,'numbertitle','off','name',text);


% 
% c=get(gcf,'children');
% % find the top of uppermost subplot
% y=0;
% for k=1:numel(c)
%     pos=get(c(k),'position');
%     y=max(y,pos(2)+pos(4));
% end
% 
%  h=annotation('textbox',[0 y 1 1-y], 'linestyle', 'none', 'string', ...
%      text, 'interpreter', 'none', 'verticalalignment', 'top', 'horizontalalignment', 'center', varargin{:});
%  
%  if nargout==0; clear h; end
%  