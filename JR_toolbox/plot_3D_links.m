function varargout = plot_3D_links(pnt, chan1, chan2, value, cfg)
% Q3D = plot_3D_links(pnt, chan1, chan2, cfg)
% - inputs
%       pnt: all channels positions (n_channels,3_xyz) (get from layout);
%       chan1: channel 1
%       chan2: channel 2
%- - output: Q3D: plot it as follow: 
%        plot3(Q3D(:,1),Q3D(:,2),Q3D(:,3)
%
% requires bezierInterp: 
% http://www.mathworks.com/matlabcentral/fileexchange/7441-bezier-interpolation-in-n-dimension-space
%
% (c) Jean-Rémi King 2011

if nargin == 4, cfg = []; end
if ~isfield(cfg,'colors'),      cfg.colors      = colormap('jet');  end     % default colormap
if ~isfield(cfg,'clim'),        cfg.clim        = [0 1];            end
if value <cfg.clim(1), value = cfg.clim(1);end
if value >cfg.clim(2), value = cfg.clim(2);end
if ~isfield(cfg,'res'),         cfg.res         = 500;              end     % arc resolutions
if ~isfield(cfg,'pltchan'),     cfg.pltchan     = 0;                end     % plot all channel
if ~isfield(cfg,'pltconnect'),  cfg.pltconnect  = 1;                end     % plot connection
if ~isfield(cfg,'loop_type'),   cfg.loop_type   = 'sphere';         end     % type of bezier
if ~isfield(cfg,'radius'),      cfg.radius      = 'sphere';         end     % type of radius (sphere, euclidian)
if ~isfield(cfg,'middle'),      cfg.middle      = 'sphere';         end     

if ~isfield(cfg,'color'),       cfg.color       = [0 0 0];          end     % line color
if ~isfield(cfg,'curve'),       cfg.curve       = 1;              end     % bezier curvature

%-- transform value into color
if size(cfg.colors,1) == 1, cfg.colors = repmat(cfg.colors,2,1); end
z = ceil((length(cfg.colors)-1)*...
    (value- cfg.clim(1))...
    /(cfg.clim(2)-cfg.clim(1)));
if z < 1, z=1;end
if z > length(cfg.colors), z = length(cfg.colors); end
color = cfg.colors(z,:);
if ~isfield(cfg,'alpha_bsl'), cfg.alpha_bsl = 1; end
if ~isfield(cfg,'alpha'), 
    cfg.alpha = cfg.alpha_bsl*(value - cfg.clim(1))/(cfg.clim(2)-cfg.clim(1)); end
if ~isfield(cfg,'linewidth_bsl'), cfg.linewidth_bsl = 10; end
if ~isfield(cfg,'linewidth'), 
    cfg.linewidth   = 1+cfg.linewidth_bsl*(value - cfg.clim(1))/(cfg.clim(2)-cfg.clim(1));
end     

 
% calculate center of head 
t=linspace(0,1,cfg.res);

% 3D Bezier Interplation
if cfg.pltchan
    scatter3(pnt(:,1),pnt(:,2),pnt(:,3));hold on;
end
%-- calculate radius
switch cfg.radius
    case 'sphere'
        if ~isfield(cfg,'Rminus'),      cfg.Rminus      = .13;              end     % R minus
        if ~isfield(cfg,'Rtimes'),      cfg.Rtimes      = 5;                end     % R times
        if ~isfield(cfg,'Rexp'),        cfg.Rexp        = 1;              end     % R exp
        R = (1+(max(sqrt(sum(pnt(chan1,:).^2)),sqrt(sum(pnt(chan2,:).^2))))*2*atan2(norm(cross(pnt(chan1,:),pnt(chan2,:))),dot(pnt(chan1,:),pnt(chan2,:)))*cfg.Rtimes)^cfg.Rexp-cfg.Rminus;
    case 'euclidian'
        if ~isfield(cfg,'Rminus'),      cfg.Rminus      = 0;              end     % R minus
        if ~isfield(cfg,'Rtimes'),      cfg.Rtimes      = 60;                end     % R times
        if ~isfield(cfg,'Rexp'),        cfg.Rexp        = 1;              end     % R exp
        R = (1+(sum((pnt(chan1,:)-pnt(chan2,:)).^2-cfg.Rminus)*cfg.Rtimes)).^cfg.Rexp-cfg.Rminus;
    otherwise
        error('unknown radius type')
end
%disp(R);
%-- calculate bezier bases
switch cfg.middle
    case 'sphere'
    [tmp tmp r1]        = cart2sph(pnt(chan1,1),pnt(chan1,2),pnt(chan1,3));     % transform into sphere coordinate
    [tmp tmp r2]        = cart2sph(pnt(chan2,1),pnt(chan2,2),pnt(chan2,3));
    [th phi tmp]        = cart2sph(mean(pnt([chan1 chan2],1)),mean(pnt([chan1 chan2],2)),mean(pnt([chan1 chan2],3))); % find middle on sphere
    [middle(1) middle(2) middle(3)] = sph2cart(th,phi,(r1+r2+r2)/3*R);        % inflate sphere
    case 'center'
        middle = mean(pnt);
    case 'origin'
        middle = [0 0 0];
    case 'manual'
        middle = cfg.middle_value;
end
control_11          = middle + (pnt(chan1,:) - pnt(chan2,:))*cfg.curve;             % create bezier controllers
control_12          = middle + (pnt(chan1,:) - pnt(chan2,:))*(1-cfg.curve);
control_21          = middle - (pnt(chan1,:) - pnt(chan2,:))*(1-cfg.curve);
control_22          = middle - (pnt(chan1,:) - pnt(chan2,:))*cfg.curve;

%-- show constructors
% scatter3(middle(1),middle(2),middle(3));
% scatter3(control_11(1),control_11(2),control_11(3), 'r');
% scatter3(control_22(1),control_22(2),control_22(3), 'g');

%-- construct bezier 1
Px=[pnt(chan1,1) control_11(1) control_12(1) middle(1)];
Py=[pnt(chan1,2) control_11(2) control_12(2) middle(2)];
Pz=[pnt(chan1,3) control_11(3) control_12(3) middle(3)];
[Q3D_1]=bezierInterp(...
    [Px(1),Py(1),Pz(1)],...
    [Px(2),Py(2),Pz(2)],...
    [Px(3),Py(3),Pz(3)],...
    [Px(4),Py(4),Pz(4)],t);
Px=[middle(1) control_21(1) control_22(1) pnt(chan2,1)];
Py=[middle(2) control_21(2) control_22(2) pnt(chan2,2)];
Pz=[middle(3) control_21(3) control_22(3) pnt(chan2,3)];
[Q3D_2]=bezierInterp(...
    [Px(1),Py(1),Pz(1)],...
    [Px(2),Py(2),Pz(2)],...
    [Px(3),Py(3),Pz(3)],...
    [Px(4),Py(4),Pz(4)],t);
Q3D = cat(1,Q3D_1(1:end-4,:),Q3D_2(4:end,:));
%-- color
if nargout == 1
    varargout = {Q3D};
else
    %-- plot connection
    if cfg.pltconnect
        h = patch(Q3D(:,1),Q3D(:,2),Q3D(:,3),'k', 'FaceAlpha', 0, 'EdgeAlpha', cfg.alpha, 'EdgeColor', color, 'LineWidth',cfg.linewidth);
%         varargout = {[Q3D_1 Q3D_2], h1, h2};
    end
end
return
