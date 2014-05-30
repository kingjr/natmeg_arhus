function [p, h, stats] = ranksum_matrix(x,y,varargin)
% [p, h, stats] = ranksum_matrix(x,y,varargin)
% same as ranksum but allows mutldimensional arrays as input
% removes NaN values in each dimension
% enable matlabpool for faster computations
% (c) JR KING, jeanremi.king+matlab@gmail.com
if size(x,2) ~= size(y,2), error('x and y should be of the same dimensions!');end

original_size = size(x);
x = reshape(x,size(x,1),[]);
y = reshape(y,size(y,1),[]);

[p h] = deal(NaN(size(x,2),1));
dim = 1;
[p(dim), h(dim) stats(dim)] = inloop(x(:,dim),y(:,dim),varargin);
parfor dim = 2:size(x,2)
    if mod(dim,size(x,2)/20)==1, fprintf('*');end
    [p(dim), h(dim) stats(dim)] = inloop(x(:,dim),y(:,dim),varargin);
end
%-- reshape data
if length(original_size) > 2
    p = reshape(p,original_size(2:end));
    h = reshape(h,original_size(2:end));
end

function [p h stats] = inloop(xt,yt,varargin)
p = 1;h = 0; stats = [];
stats.zval = [];
stats.ranksum= [];

%-- deal with nan
xt(isnan(xt)) = [];
yt(isnan(yt)) = [];
if ~isempty(xt) && ~isempty(yt)
    if ~isempty(varargin)
        [p h stats] = main(xt,yt,varargin);
    else
        [p h stats] = main(xt,yt);
    end
end
return


function [p, h, stats] = main(x,y,varargin)
add = [];
try
    if ~isempty(varargin)
        for arg = 1:length(varargin)
            add = [add ', ' varargin{arg}];
        end
    end
    eval(['[p h stats] = ranksum(x,y' add ');']);
catch
    eval(['[p h stats] = ranksum(x,y);']);
end
return