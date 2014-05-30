function data = hatchUnsig2(x,y,data,hc,cfg)
%  hatchUnsig2(x,y,data,hc,cfg)
% function imagesc data and add hatches where h == 0
% x and y are used as the ticks labels
% can do it repeatively if data and h are structures
% -- cfg.
%       linesDensity
%       color
%       linewidth

%-- deal with inputs
if nargin <= 3,data = x;h = y; x = 1:size(data,1);y = 1:size(data,2);end
if nargin == 2 || nargin == 4, cfg.tmp = [];end
if ~isfield(cfg,'linesDensity'),    cfg.linesDensity    = 500;  end         % hatches density
if ~isfield(cfg,'linewidth'),       cfg.linewidth       = 1;    end         % hatches width
if ~isfield(cfg,'color'),           cfg.color           = 'k';  end         % hatches color
if not(iscell(hc)), hc = {hc};                                  end         % for multiple hatches layers only

%-- if labels don't correspond to matrix
if length(x) ~= size(data,2) || length(y) ~= size(data,1)
    disp('error label size')
    x = x(1:size(data,2));
    y = y(1:size(data,1));
end

%-- deal with borders
x(end+1)        = x(end) + (x(end) - x(end-1));
y(end+1)        = y(end) + (y(end) - y(end-1));
x(1:end-1)      = x(1:end-1)-(x(2:end) - x(1:(end-1)))./2;
y(1:end-1)      = y(1:end-1)-(y(2:end) - y(1:(end-1)))./2;
data(end+1,:)   = data(end,:);
data(:,end+1)   = data(:,end);
for c = 1:length(hc), hc{c}(end+1,:) = hc{c}(end,:);hc{c}(:,end+1) = hc{c}(:,end);end

hold on;
surf(x,y,ones(size(data,1),size(data,2))-.5,data,'edgealpha',0);   % surf all image

%-- build up hatches: create lines
h_x   = [-2:size(data,1)/cfg.linesDensity:1;-1:size(data,1)/cfg.linesDensity:2];
h_y   = [zeros(1,length(h_x))-1; ones(1,length(h_x))+2];
h_z   = ones(2,length(h_x));
if ~iscell(hc)
    hc = {hc};
end

for c = 1:length(hc)
    
    h = hc{c};
    %-- removes non hatched values
    data2           = data;
    data2(h == 0)   = NaN;
    
    %-- turn hatches
    [th r]    = cart2pol(h_x,h_y);
    [h_x2 h_y2] = pol2cart(th + pi*(c-1)/length(hc),r);
    
    %-- scale hatches
    h_x2 = (max(x)-min(x)).*h_x2;
    h_y2   = (max(y)-min(y)) .* h_y2;
    
    %-- print figure
    plot3(h_x2,h_y2,h_z+c-1,cfg.color,'linewidth',cfg.linewidth);                  % plot lines
    surf(x,y,ones(size(data,1),size(data,2))+c-.5,data2,'edgealpha',0);  % surf non hatched image
    axis([min(x) max(x) min(y) max(y)-.5 0 c+1]);                        % set axis
    
end
return