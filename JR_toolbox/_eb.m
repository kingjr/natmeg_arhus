function eb(x,y,options)

%function eb(x,y,options)
%------------------------------------------------------------------
% plot error y as function of x, with a transparent errorbar
% options is optional, and defines, in a structure
%   options.color       
%   options.alpha
%   options.edgeAlpha
%   options.lineAlpha
%   options.lineWidth
%   options.eb:        'std' or 'eb'
% (c) JeanRémi King, all rights reserved
%------------------------------------------------------------------
if nargin == 2,options.tmp = [];end
options.color     = getoption(options,'color', [0 0 1]);
options.alpha     = getoption(options,'alpha', .5);
options.edgeAlpha = getoption(options,'edgeAlpha', 0);
options.lineAlpha = getoption(options,'lineAlpha', 1);
options.lineWidth = getoption(options,'lineWidth', 1);
options.eb        = getoption(options,'eb', 'eb');

%-- error bar
switch options.eb
    case 'std'
        e = std(y);
    case 'eb'
        e = std(y) ./ sqrt(size(y,1));
end

hold on;
% plot fill vector
linehandle=fill(...
    [x fliplr(x)],...
    [mean(y) fliplr(mean(y))],...
    options.color);
set(linehandle,...
    'lineWidth',options.lineWidth,...
    'EdgeColor',options.color,...
    'FaceAlpha',options.alpha,...
    'EdgeAlpha',options.lineAlpha...
    );

fillhandle=fill(...
    [x fliplr(x)],...
    [mean(y)-e fliplr(mean(y)+e)],...
    options.color);

% set transparency
set(fillhandle,...
    'EdgeColor',options.color,...
    'FaceAlpha',options.alpha,...
    'EdgeAlpha',options.edgeAlpha...
    );

hold off;
end