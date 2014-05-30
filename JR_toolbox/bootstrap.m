function [h p] = bootstrap(x,m,b,alpha,n,void)
% FUNCTION BOOSTRAP
% bootstrap(x,m,b,alpha,n,parametric,void)
% do a resampling of x, according to axis 1, n times, evaluate mean, and
% check whether it's significantly higher than m
% -- input
%      x:           matrix where lines are different instances, and columns
%                   dimensions
%      m:           mean to be compared to
%      b:           (optional) number of resampling
%      alpha:       (optional) alpha value for test
%      n:           (optional) bonferroni correction
%      parametric:  (optional) h = h_non-parametric * h_parametric
%      void:        (optional) show progressbar
% -- output
%      h:           hypothesis verified up to alpha
%      p:           p value
%
% parametric depends on stats toolbox
% void necessitate progressbar.m
%--------------------------------------------------------------------------
% (c) JeanRémi KING: jeanremi.king+matlab@gmail.com, all rights reserved
%--------------------------------------------------------------------------
if nargin <= 2
    b = 1000; % number of boostrap loop
end
if nargin <= 3
    alpha = .05; % alpha value
end
if nargin <= 4
    n = numel(x) / size(x,1) ; % bonferronni correction
end
if nargin <= 5
    parametric = 0;
end
if nargin <= 6
    void = 0;
end
alpha = alpha / n;

% reshape x into 2D matrix
xsize = size(x);
x = reshape(x,size(x,1),[]);


d = NaN .* zeros(b,size(x,2));
for ii = 1:b % boostrap n times
    if void == 1
        try progressbar(ii,b);% need progressbar
        catch exception
        end
    end
    % make new distribution
    d(ii,:) = mean(x(randi(size(x,1),1,size(x,1)),:));
end
% t-test whether it's higher than m
p = (sum(d - m < 0) / b);
h = p < alpha;

%return array as it was
try
    p = reshape(p,xsize(2:end)) .* reshape(~isnan(squeeze(mean(d))),xsize(2:end));
    h = reshape(h,xsize(2:end)) .* reshape(~isnan(squeeze(mean(d))),xsize(2:end));
end

return