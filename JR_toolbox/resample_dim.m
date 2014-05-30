function X = resample_dim(X,P,Q,dim,memory)
% Y = resample_dim(X,P,Q,dim);
% resample n-dimensional matrix function
% (c) JeanRémi King 2012
if nargin == 4, memory = 'high'; end
dims    = size(X);
X       = permute(X,[dim setdiff(1:length(dims),dim)]);
index   = [dim setdiff(1:length(dims),dim)];
X       = reshape(X,size(X,1),[]);
switch memory
    case 'high'
        X       = resample(X,P,Q);
    case 'low'
        for d = size(X,2):-1:1 % to automatically initialize at largest size
            Xr(:,d) = resample(X(:,d),P,Q);
        end
        X = Xr; clear Xr;
end
X       = reshape(X,[size(X,1), dims(setdiff(1:length(dims),dim))]);
[unused index_sorted] = sort(index);
X       = permute(X, index_sorted);
