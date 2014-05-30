function my_plot_topo_fast(w,cfg)
% my_plot_topo_fast(values,cfg)
% fast topographical plot based on scatter rather than interpolated images.
% (c) JeanRémi King 2012
% jeanremi.king [at] gmail.com
if nargin == 1, cfg = []; end
if ~isfield(cfg,'map'),     cfg.map     = colormap('jet');  end             % color code
if ~isfield(cfg,'size'),    cfg.size    = 20;               end             % scatter size
if ~isfield(cfg,'chans'),   cfg.chans   = 1:length(w);      end             % selected channels in layout
if ~isfield(cfg,'layout'),                                                  % prepare layout
    switch length(w) 
        case 306,   cfg.layout  = 'neuromag306all.lay'; cfg.layout = ft_prepare_layout(cfg);
        case 206,   cfg.layout  = 'neuromag306grad.lay';cfg.layout = ft_prepare_layout(cfg);
        case 204,   cfg.layout  = 'neuromag306grad.lay';cfg.layout = ft_prepare_layout(cfg);
        case 104,   cfg.layout  = 'neuromag306mag.lay'; cfg.layout = ft_prepare_layout(cfg);
        case 102,   cfg.layout  = 'neuromag306mag.lay'; cfg.layout = ft_prepare_layout(cfg);
        case 256,   cfg.layout  = load('my_EGI_net.mat', 'layout'); cfg.layout=cfg.layout.layout; 
        otherwise
            error('specify cfg.layout!');
    end
end
if ~isfield(cfg,'colors')
if ~isfield(cfg,'clim'), cfg.clim = [min(w) max(w)]; end
%-- normalize color according to values
w = ceil((w-cfg.clim(1))/(cfg.clim(2)-cfg.clim(1)).*size(cfg.map,1));
w(w<1) = 1;
w(w>size(cfg.map,1)) = size(cfg.map,1);
cfg.colors = cfg.map(w,:);
end

%-- plot
scatter(cfg.layout.pos(cfg.chans,1),cfg.layout.pos(cfg.chans,2),cfg.size,cfg.colors,'filled');
axis off image;