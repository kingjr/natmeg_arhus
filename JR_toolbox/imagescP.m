function imagescP(X,Y,data,testsNb,alpha)

if nargin == 1
    data = X; 
elseif nargin <= 3
    testsNb = 1;
    alpha = .01;
end
load('MyColormaps.mat');
alpha = alpha / testsNb;%alpha bonferroni corrected
meandata = squeeze(mean(data));%mean
s=squeeze(nanstd(data) / size(data,1));%standard error
% MSE standard error of the mean: 2 tailed margins:
MSE = norminv(1 - alpha / 2);
% MSE map show where it's not significant:
MSE_map = ((((meandata + MSE .* s) >= .5) .* (meandata <= .5)) + ...
    (((meandata - MSE .* s) <= .5) .* (meandata >= .5)));

% imagesc(X,Y,meandata - MSE_map);
% axis([xlim ylim zlim -1 1]);
% set(gcf,'colormap', 'mycolormap');

imagesc(X,Y,meandata .* ~MSE_map);
axis([xlim ylim zlim 0.4 .75]);
set(gcf,'colormap', myhot);

return