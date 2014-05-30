function [CData handle]= topo_plot_color(X,cfg)
% CData = topo_plot_color(X,cfg);
handle = [];
if ~isfield(cfg, 'plot'), cfg.plot = 1;end
if ~isfield(cfg, 'rotate_color'), cfg.rotate_color = 1;end
% build templates
for group = 1:size(X,1)
    if group == 1 || nansum(abs(X(group,:))) > 0
        figure();clf;set(gcf,'color','w');
        my_plot_topo(255*X(group,:)',cfg);
        axis([-.5 .5 -.5 .5 0 255 0 255]);colormap('gray');
        c = get(gca, 'children');%extract interpolation image
        if group == 1
            C = get_CData(c);
            CData = repmat(C,[1 1 3]);
            CData(~isnan(CData)) = 0;
            CData(:,:,group) = C;
        else
            CData(:,:,group) = get_CData(c);
        end
        close;
    end
end
%-- make sure correct lims
CData(CData<0) = 0;
CData(CData>255) = 255;
%
% if cfg.rotate_color
%
%     real = CData;
%     real(~isnan(real)) = 1;
%     CDataRGB = CData;
%     CDataRGB(isnan(CData)) = 0;
%     CDataRGB = rotate_color(CDataRGB, -1/6);
%     imagesc(CDataRGB./255);
%     CData = 255*CDataRGB .* real;
% end

% plot
if cfg.plot
    my_plot_topo(zeros(256,1),cfg);
    c = get(gca, 'children');%extract interpolation image
    c = set_CData(c,CData./255);
    %     set(c(2), 'visible', 'off');
    axis image off;
    handle = gca;
end

end

function C = get_CData(c)
C = [];
for f = 1:length(c)
    properties = get(c(f));
    if isfield(properties,'CData') && ~isempty(properties.CData)
        C = properties.CData;
        return
    end
end
end


function c = set_CData(c,C)
for f = 1:length(c)
    properties = get(c(f));
    if isfield(properties,'CData') && ~isempty(properties.CData)
        set(c(f), 'CData', C, 'CDataMapping', 'direct');
        return
    end
end
end
