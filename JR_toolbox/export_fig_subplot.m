function export_fig_subplot(name,varargin)
if nargin == 1, varargin = {}; end
for ii = 1:2:length(varargin);
    eval([varargin{ii} '=varargin{ii+1};']);
end
if ~exist('fig', 'var'), fig = gcf; end
if ~exist('ext', 'var'), ext = 'png'; end
children = get(fig, 'children');
for c = 1:length(children)
    f = figure(9999);clf;set(gcf,'color','w');
    copyobj(children(c),figure(9999));
    set(gca,'position',[.05 .05 .9 .9]);
    export_fig([name '_' num2str(c) '.' ext]);
end