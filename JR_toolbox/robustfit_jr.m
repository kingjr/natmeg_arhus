function [brob, rob_stats] = robustfit_jr(x,y,varargin)
% [brob, rob_stats] = robustfit_jr(x,y,varargin)
% apply robustfit and also outputs R²
% see robustfit help
% script adapted from:
% http://www.mathworks.fr/support/solutions/en/data/1-CMABGO/index.html?product=ST&solution=1-CMABGO
if nargin == 2, varargin = {}; end
[bls,~,~,~,linreg_stats] = regress(y,[ones(size(x)) x]);
[brob, rob_stats] = robustfit(x,y,varargin{:});

if 0
scatter(x,y,'filled'); grid on; hold on
plot(x,bls(1)+bls(2)*x,'r','LineWidth',2);
plot(x,brob(1)+brob(2)*x,'g','LineWidth',2)
legend('Data','Ordinary Least Squares','Robust Regression')
end
rob_stats.rquare_linreg = linreg_stats(1);
rob_stats.rsquare_robustfit = corr(y,brob(1)+brob(2)*x)^2;