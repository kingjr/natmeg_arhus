function [curve auc ]=myroc(distr1,distr2)
%[curve auc stats]=myroc(distr1,distr2,[cfg])
% calculates, from the distributions (on 1 dimension) of two populations
% the discrete ROC curves and the Area Under the Curve
% (c) JeanRémi King 2011
%%
both    = sort(cat(1,distr1,distr2));                                       % calculate curve
curve   = [0,0];                                                            % starting point
if length(both) < 1000, 
    x = 1:length(both);
else
    x = linspace(1,length(both),1000);
end
ii = 0;
for point = 1:length(x)
    ii = ii+1;
    curve(ii,1)  = length(find(distr1*10^15<=(both(round(x(point)))*10^15)))/size(distr1,1);   % probability of distr1 under a given point
    curve(ii,2)  = length(find(distr2*10^15<=(both(round(x(point)))*10^15)))/size(distr2,1);
end
auc = trapz(curve(:,1),curve(:,2));                                         % calculate area under the curve
return