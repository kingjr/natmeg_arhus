function  [handle_line handle_fill] = plot_eb(X,Y,cfg)
% plot_eb(X,Y,[cfg])
% plot error bar
% (c) JeanRémi King, jeanremi.king+matlab[at]gmail.com

if nargin < 3
    cfg = [];
end
evalc('if length(cfg) == 3, cfg.color = cfg; end');
if ~isfield(cfg, 'color'),  cfg.color   = [.2 .3 .9]; end
if ~isfield(cfg, 'type'),   cfg.type    = 'sem'; end
if ~isfield(cfg, 'rmnan'),  cfg.rmnan   = true; end
if ~isfield(cfg, 'FaceAlpha'),  cfg.FaceAlpha   = .5; end
if ~isfield(cfg, 'EdgeAlpha'),  cfg.EdgeAlpha   = .5; end


if cfg.rmnan
    X = X(~isnan(nanmean(Y)));
    Y = Y(:,~isnan(nanmean(Y)));
end

m=nanmean(Y);

index_1=1:length(m);
index_2=index_1(end:-1:1);
% s = nanstd(Y) / (size(Y,1) % does not take NaN nto account


s = NaN;
switch cfg.type
    case 'sem'
        for x = 1:length(m)
            s(x)=nanstd(Y(:,x)) ./ sqrt(size(Y,1) - sum(isnan(Y(:,x))));
        end
        area = [m-s m(index_2)+s(index_2)];
    case 'std'
        for x = 1:length(m)
            s(x)=nanstd(Y(:,x));
        end
        area = [m-s m(index_2)+s(index_2)];
    case 'prctile'
        for x = 1:length(m)
            s1(x)=prctile(Y(:,x),10);
            s2(x)=prctile(Y(:,x),90);
        end
        area = [s1 s2(index_2)];
end


hold on;
handle_line=plot(X,m,'linewidth',2,'Color',cfg.color);
try
    handle_fill=fill([X flipdim(X,2)],area,cfg.color);
    set(handle_fill,'edgecolor',cfg.color, 'FaceAlpha', cfg.FaceAlpha, 'EdgeAlpha', cfg.EdgeAlpha);
end

end


