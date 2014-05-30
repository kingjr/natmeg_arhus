function r = getRank(X)
[Y I] = sort(X);
r(I) = 1:length(X);