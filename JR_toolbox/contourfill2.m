function [poly h] = contourfill2(xs,ys,data,level, varargin)
hold on;
m = max(data(:));
for ii = length(level):-1:1
    [poly{ii} h{ii}] = contourfill(xs,ys,data,[level(ii) m], varargin{:});
end