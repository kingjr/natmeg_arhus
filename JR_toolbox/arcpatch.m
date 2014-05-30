function [varargout] = arcpatch(x,y,r,phi)

% ARCPATCH creates a semicircle wedge as a MATLAB patch. Useful for
%   highlighting the angle between two lines, low level pie chart
%   creation, etc.
% 
%   USES:
%       arcpatch(x,y,r,phi)
%       arcpatch(x,y,r,[phi0 phi1])
%       h = archpatch(...)
% 
% 
%   arcpatch(x,y,r,phi)             creates an arcpatch with center (x,y),
%                                   radius r, and from zero to phi degrees
%   arcpatch(x,y,r,[phi0 phi1])     creates an arcpatch from phi0 to phi1
%                                   degrees
%   h = arcpatch(...)               returns the handle of the patch object
% 
%
%   See also PATCH (built-in)
%   See also CIRCULARC (File ID#4082) available from the File Exchange
%
% 
%   Author:
%   Todd C Pataky (0todd0@gmail.com)  ['zero' todd 'zero'@gmail.com]
%   01-October-2007


% EXAMPLES:
% (1) Basic use:
%       h = arcpatch(0,0,1,60);
%       set(h,'facecolor','c','edgecolor','b')
%       axis equal
%
% (2) Indicating the angle between two lines:
%       lh(1) = line([0 1],[0 0]);
%       lh(2) = line([0 1],[0 1]);
%       set(lh,'color','k','linewidth',2)
%       h = arcpatch(0,0,0.4,45);
%       set(h,'facecolor',0.5*[1 1 1],'edgecolor','k')
%       th = text(0.2,0.1,'\theta');
%       set(th,'color','k','fontsize',20)
%       axis equal off
%
% (3) Low level pie chart creation:
%       h = zeros(1,5);
%       h(1) = arcpatch(0,0,1,30);
%       h(2) = arcpatch(0,0,1,[30 120]);
%       h(3) = arcpatch(0,0,1.1,[120 180]);
%       h(4) = arcpatch(-0.1,0,0.9,[180 210]);
%       h(5) = arcpatch(0,0,1.2,[200 370]);
%       set(h(1),'facecolor','k')
%       set(h(2),'facecolor','c')
%       set(h(3),'facecolor','b')
%       set(h(4),'facecolor','m')
%       set(h(5),'facecolor','r')
%       set(h,'edgecolor','none','facealpha',0.5)
%       axis equal off
%



%(0) Identify which of the two USES (above) has been specified
switch length(phi)
    case 1      %first of the USES
        [phi0,phi1] = deal(0,phi);   
    case 2      %second of the USES
        [phi0,phi1] = deal(phi(1),phi(2));
    otherwise
        error('Argument ''phi'' must have either one or two elements.')
end
clear phi



%(1) Define angles at which patch vertices will appear:
step = 2;               %one vertex for every 2 degrees; adjust this parameter for smoother/coarser edges
phi = phi0:step:phi1;   %angles (degrees) at which patch vertices will appear
if phi(end)~=phi1
    phi = [phi  phi1];  %append phi1 if 'step' is not a factor of (phi1-phi0) 
end
phi = pi/180*phi';      %angles (radians)


%(2) Define patch vertices and faces
vert = zeros(length(phi)+1,3);                      %extra vertex for the arc center
vert(1,:) = [ x y 0 ];                              %arc center
vert(2:end,1:2) = [x+r*cos(phi)   y+r*sin(phi)];    %arc radius
fac = [1:size(vert,1) 1];                           %rejoin the arc radius with the center


%(3) Create patch
h = patch('vertices',vert,'faces',fac);


if nargout==1
    varargout{1} = h;  %return patch handle
end

    

