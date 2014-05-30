function [p h stats] = signrank_matrix(x,y,varargin)

s = size(x);

x = reshape(x,s(1),[]);
y = reshape(y,s(1),[]);
for dim = size(x,2):-1:1
    [p(dim) h(dim) stats(dim)] = signrank(x(:,dim),y(:,dim),varargin{:});
end
if length(s)>2
    p = reshape(p,s(2:end));
    h = reshape(h,s(2:end));
    stats = reshape(h,stats(2:end));
end