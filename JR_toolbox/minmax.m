function var = minmax(x,dim)
% var = minmax(x,dim)
% jeanremi [dot] king [at] gmail [dot] com
if nargin==1,
    if sum(size(x))>2
        dim=find(size(x)>1,1);
    else
        dim = 2;
    end
end

var = cat(dim,min(x,[],dim),max(x,[],dim));