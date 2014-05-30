function Y = propGraph(input,options)
% propGraph(input,res,ths)
% plot distribution of input(:,t) < [res]
% -- input: 2D matrix, sum across first dimention
% -- options.res: resolution of the plot
% -- options.ths: plot specific line thresholds
% -- options.skip: don't calculate proportion
% -- outputs matrix Y

if nargin == 1
    options.tmp = [];
end
res = getoption(options,'res',100);
ths = getoption(options,'ths',[]);
skip = getoption(options,'skip',0);
c = getoption(options,'colormap', colormap);
if not(skip)
    Y(1,:) =  squeeze(sum((input < 1/res)))./size(input,1);
    lineCol{1} = c(round(1/res*size(colormap,1)),:);
    for r = 2:res
        Y(r,:) = squeeze(sum((input <= r/res) .* (input > (r-1)/res)))./size(input,1);
        
    end
else
    Y = input;
end
h = area(Y');
hold on;
for r = 1:res
    if ismember(r/res,ths)
        plot(sum(Y(1:r,:),1),'color',[0 0 0],'linewidth', 3);
    else
        plot(sum(Y(1:r,:),1),'color',c(round(r/res*size(colormap,1)),:));
    end
end
axis([xlim 0 1]);
hold off;
