function handle = ft_topoplotPlanars(cfg,data,varargin)

if ~isfield(cfg, 'layout'), error('need cfg.layout');end
if ischar(cfg.layout), cfg.layout = ft_prepare_layout(cfg); end

% select data
selchan = @(c,s) find(cell2mat(cellfun(@(x) ~isempty(strfind(x,s)),c,'uniformoutput', false))==1);
try
channels = intersect(...
    selchan(cfg.layout.label,'+'),...
    selchan(cfg.layout.label,'MEG'));
catch e
    warning('you have to combine gradiometers first!')
    return
end

% time average
avg = mean(data.avg(channels,find(data.time>=cfg.xlim(1),1):find(data.time>=cfg.xlim(2),1)),2);

% scale
if ~isfield(cfg, 'zlim'), cfg.zlim = 1./[min(avg(:)) max(avg(:))]; end
[th r] = cart2pol(real(avg),imag(avg));
r(r(:)<cfg.zlim(1)) = cfg.zlim(1); % min plot value
r(r(:)>cfg.zlim(2)) = cfg.zlim(2); % max plot value
r = (r-cfg.zlim(1))./(cfg.zlim(2)-cfg.zlim(1)); % value between 0 and 1
[avg(:,1) avg(:,2)] = pol2cart(th,r);
avg = avg(:,1)+1i*avg(:,2);

% plot type
if ~isfield(cfg,'plot_type'), cfg.plot_type = 'quiver_angle';end
% arrow_size
if ~isfield(cfg,'size'), cfg.size = .1; end

holdon = ishold; hold on;
switch cfg.plot_type
    case 'quiver'
        % color
        if ~isfield(cfg,'color'), cfg.color = [0 0 0]; end
        % quiver
        hold on;
        for c = 102:-1:1
            handle{c}= quiver(...
            cfg.layout.pos(channels(c),1),cfg.layout.pos(channels(c),2),...
            real(avg(channels(c))), imag(avg(channels(c))),...
            cfg.size,'color',cfg.color, 'linewidth',2);
        end
        % axes
        axis([-1 1 -1 1] * .6);axis square off;
    case 'quiver_norm'
        % normalize data into unit vectors
        [x,y] = pol2cart(angle(avg),.5);
        xy = x+1i*y;
        if ~isfield(cfg,'color'), cfg.color = colormap('jet'); end
        for c = 102:-1:1
            handle{c}=quiver(...
                cfg.layout.pos(channels(c),1),cfg.layout.pos(channels(c),2),...
                real(xy(c)), imag(xy(c)),...
                cfg.size,...
                'color', cfg.color(round(abs(avg(channels(c)))*(size(cfg.color,1)-1))+1,:),...
                'linewidth',2);
        end
        axis([-1 1 -1 1] * .6);axis square off;
    case 'quiver_angle'
        holdon = ishold; hold on;
        % compute hsv color map
        hue         = (pi+angle(avg))/(2*pi);
        saturation  = abs(avg);
        hsv         = hsv2rgb(cat(2,hue,saturation,ones(length(avg),1)));
        % plot
        for c = 102:-1:1
            handle{c}=quiver(...
                cfg.layout.pos(channels(c),1),cfg.layout.pos(channels(c),2),...
                real(avg(channels(c))), imag(avg(channels(c))),...
                cfg.size,...
                'color', hsv(c,:),...
                'linewidth',2);
        end
        axis([-1 1 -1 1] * .6);axis square off;
    case 'interp'
        % interpolation
        if ~isfield(cfg,'resolution'), cfg.resolution = 100; end
        % define interp
        xi = linspace(...
            min(cfg.layout.pos(channels,1)),...
            max(cfg.layout.pos(channels,1)),...
            cfg.resolution);
        yi = linspace(...
            min(cfg.layout.pos(channels,2)),...
            max(cfg.layout.pos(channels,2)),...
            cfg.resolution);
        [xi yi] = meshgrid(xi,yi);
        interpolate = @(z) griddata(...
            cfg.layout.pos(channels,1),cfg.layout.pos(channels,2),...
            z,...
            xi,yi, 'cubic');
        % transform data into RGB color
        hue         = (pi+angle(avg))/(2*pi);
        saturation  = abs(avg(1:102));
        rgb         = hsv2rgb(cat(2,hue,saturation,ones(length(avg),1)));
        
        z = cat(3,...
            interpolate(rgb(:,1)),...
            interpolate(rgb(:,2)),...
            interpolate(rgb(:,3)));
        z(z(:)<=0) = eps; % correct matlab bug
        z(z(:)>=1) = 1-eps;
        
        imagesc(z);
        axis image off;set(gca,'ydir','normal');
end
if ~holdon, hold off; end
