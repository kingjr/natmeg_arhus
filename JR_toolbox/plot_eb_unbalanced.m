function plot_eb_unbalanced(Xr,Yr,cfg)
% idem as myeb(X,Y,cfg.color) but for unebalanced data

if nargin < 3 ,    cfg = []; end
if ~isfield(cfg, 'color'),  cfg.color   = [.2 .3 .9];     end % colors
if ~isfield(cfg, 'bin'),    cfg.bin     = size(Xr,2)/10;   end % number of bins

%-- cuts into bin
X = (min(Xr(:))-.5*(max(Xr(:))/cfg.bin)):(max(Xr(:))/cfg.bin):(max(Xr(:)));
if X(end) ~= max(Xr(:)), X(end+1) = max(Xr(:)); end

for bin = 1:(length(X)-1)
    if bin == (length(X)-1)
        ind = (Xr > X(bin)).* (Xr < X(bin+1));
    else
        ind = (Xr > X(bin)).* (Xr < X(bin+1));
    end
    for ii = 1:size(Yr,1)
        Y(ii,bin) = nanmean(Yr(ii,find(ind(ii,:))));
    end
end
%     
X = X(2:end) - (X(2)-X(1))/2;


m=nanmean(Y);
s=nanstd(Y) / size(Y,1); 
index_1=1:length(m);
index_2=index_1(end:-1:1);

hold on; 
g=fill([X flipdim(X,2)],[m-s m(index_2)+s(index_2)],cfg.color);
plot(X,m,'linewidth',2,'Color',cfg.color);
set(g,'edgecolor',cfg.color, 'FaceAlpha', 0.5, 'EdgeAlpha', 0); 
