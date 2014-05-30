function colors = mycolormap(n,name)
%colors = mycolormap(n,name)
if nargin == 1,name = 'jet';end
colors = colormap(name);
colors = colors(round(linspace(1,size(colors,1),n)),:);
end