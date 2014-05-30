function varargout=plot_arc(P1,P2,R, cfg, varargin)

if  nargin <= 5, cfg = [];                  end
if  nargin <= 6, options = {};
else             options = varargin(6:end); end
% parameters
if ~isfield(cfg, 'res'), cfg.res = 60; end % circle resolution
if ~isfield(cfg, 'arc'), cfg.arc = pi; end % circle perimeter angle

% find center of ellipse
xyz = mean([P1 ; P2]);

% find first radius
r(1) = sqrt(sum((P2-P1).^2))/2;
r(2) = R; % radius 2

% find angle
[th tmp1 tmp2] = cart2pol(P1(1)-P2(1), P1(3)-P2(3), P1(2)-P2(2));
thetaphi(1) = 180/pi*-th;
[th tmp1 tmp2] = cart2pol(P1(3)-P2(3), P1(2)-P2(2), P1(1)-P2(1));
psi = 180/pi*-th;

thetaphi(2) = 0;

% uses 60 intervals
t = linspace(0, cfg.arc, cfg.res)';

% polyline approximation of ellipse, centered and parallel to main axes
x       = r(1) * cos(t);
y       = r(2) * sin(t);
z       = zeros(length(t), 1);
base    = [x y z];

% compute transformation from local basis to world basis
trans   = localToGlobal3d(xyz(1), xyz(2), xyz(3), thetaphi(1), thetaphi(2), psi);

% transform points composing the ellipse
ellipse = transformPoint3d(base, trans);

% draw the curve
h = drawPolyline3d(ellipse, options{:});


if nargout > 0
    varargout = {h};
end
