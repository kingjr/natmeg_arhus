function handles = myboxplot(x,y,varargin)

if nargin == 2
    varargin = {};
end
groups = unique(y);

index = (find(ismember(varargin(1:2:end),'widths'))-1)*2+1;
if ~isempty(index), widths = varargin{index+1};else widths = .5;end
h = ishold;
hold on;    
for group = length(groups):-1:1
    height = prctile(x(y==groups(group)),[25 75]);
    handles_fill(group) = fill(group + [-1 1 1 -1]*widths/2,height([1 1 2 2]),[0 0 0],'edgecolor','none');
end

handles = boxplot(x,y,varargin{:});
if h
    hold on;
else
    hold off;
end
handles = cat(1,handles,handles_fill);