function handle = scattereb(x,y,color,linewidth)
% function handle = scattereb(x,y,<color>,<linewidth>)
% plots mean and SEM of y a function of x (SEM on both axes)
% X and Y are two 1-D vectors, where lines are entities
% depends on stats toolbox
% (c) Jean-Remi KING, all rights reserved, jeanremi.king+matlab@gmail.com
if nargin < 3
    color = 'k';
end
if nargin < 4
    linewidth = 1;
end

x(isnan(x)) = [];
y(isnan(y)) = [];

x_m = mean(x);
y_m = mean(y);
x_eb = std(x) / sqrt(length(x));
y_eb = std(y) / sqrt(length(y));
if isstr(color)
    handle = plot([x_m - x_eb, x_m + x_eb; x_m, x_m]', [y_m, y_m; y_m - y_eb, y_m + y_eb]', color, 'linewidth', linewidth);
elseif isvector(color)
    handle = plot([x_m - x_eb, x_m + x_eb; x_m, x_m]', [y_m, y_m; y_m - y_eb, y_m + y_eb]', 'color',color, 'linewidth', linewidth);
end
return