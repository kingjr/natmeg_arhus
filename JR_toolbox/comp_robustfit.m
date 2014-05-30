function [t p df] = comp_robustfit(XY1,XY2)

%-- get robust regression parameters
[B1,STATS1]   = ROBUSTFIT(XY1(:,1),XY1(:,2));
[B2,STATS2]   = ROBUSTFIT(XY2(:,1),XY2(:,2));
%-- standard error of the difference
S = sqrt(STATS1.se(2)^2 + STATS2.se(2)^2);
%-- ttest
t = (B1(2) - B2(2))/S;
%-- degrees of freedom
df= size(XY1,1)+size(XY2,1)-4;
%-- pvalue
p = tcdf(t,df);
return