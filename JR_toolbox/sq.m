function x = sq(x)
% fast and efficient squeeze
x = squeeze(x);
while size(x,1) == 1 && numel(x)>1
    x = squeeze(reshape(x,[size(x,2:length(size(x))),1]));
end