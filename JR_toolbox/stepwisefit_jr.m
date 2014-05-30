function [b,se,pval,inmodel,stats,n,h] = stepwisefit_jr(allx,y,varargin)
% [b,se,pval,inmodel,stats,n,h] = stepwisefit_jr(x,y,varargin)
% applies stepwisefit (see help) and also return adjusted R-square
% % from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/112926
if nargin == 2, varargin ={};end
[b,se,pval,inmodel,stats,n,h]=stepwisefit(allx,y, varargin{:});
stats.R2 = 1 -(stats.SSresid/stats.SStotal)*((stats.dfe+stats.df0)/(stats.dfe));