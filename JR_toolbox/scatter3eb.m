function scatter3eb(x,y,z,color,linewidth)
% function scattereb(x,y,z,<color>,<linewidth>)
% plots mean and SEM of y a function of x (SEM on both axes)
% X and Y are two 1-D vectors, where lines are entities
% depends on stats toolbox
% (c) Jean-Remi KING, all rights reserved, jeanremi.king+matlab@gmail.com
if nargin < 4
    color = 'k';
end
if nargin < 5
    linewidth = 1;
end

x(isnan(x)) = [];
y(isnan(y)) = [];
z(isnan(z)) = [];

x_m = mean(x);
y_m = mean(y);
z_m = mean(z);
x_eb = std(x);
y_eb = std(y);
z_eb = std(z);
if 0
    x_eb = x_eb/ sqrt(length(x));
    y_eb = y_eb / sqrt(length(y));
    z_eb = z_eb / sqrt(length(z));
end
points = [xm-eb,ym-eb,zm-eb];
if isstr(color)
    
    plot3(...
        [x_m - x_eb, x_m + x_eb; x_m, x_m]', ...
        [y_m, y_m; y_m - y_eb, y_m + y_eb]', ...
        [z_m, z_m; z_m - z_eb, z_m + z_eb]', ...
        color, 'linewidth', linewidth);
elseif isvector(color)
    plot3(...
        [x_m - x_eb, x_m + x_eb; x_m, x_m]',...
        [y_m, y_m; y_m - y_eb, y_m + y_eb]', ...
        [z_m, z_m; z_m - z_eb, z_m + z_eb]', ...
        'color',color, 'linewidth', linewidth);
end
return