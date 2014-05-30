function cfg = myvenn(data,cfg)
% plot 2 venn circles and in an overall population
% cfg = myvenn(data,cfg)
% necessitates circle.m to draw circle
% (c) JeanRémi King, jeanremi.king+matlab@gmail.com
% adapted from venn.m 
if ~isfield(cfg, 'resolution'), cfg.resolution  = 50;               end
if ~isfield(cfg, 'style'),      cfg.style       = {'k-';'r-';'b-'}; end
if ~isfield(cfg, 'position'),   cfg.position    = [0 0];            end

%-- transform data for venn2 format

%-- radius
r1 = sqrt( (data(1)+data(2)+data(3)+data(4))/pi );
r2 = sqrt( (data(2)+data(3)+data(4)+data(4))/pi );
r3 = sqrt( (data(3)+data(3)+data(4)+data(4))/pi );

%-- geometrical solution of venn problem for 3 circles
y = ( dist_A_C^2 - dist_B_C^2 + dist_A_B^2 ) / 2 / dist_A_B;
size_x = max( r1 + dist_A_B + r2, 2*r3 );
size_y = max( r1, r2 ) + sqrt( dist_A_C^2 - y^2 ) + r3;

%-- find the circle centers
center1 = [r1,              max( r1, r2 )];
center2 = [r1 + dist_A_B    center1(2)];
center3 = [r1 + y,          center1(2) + sqrt( dist_A_C^2 - y^2 )];

%draw the circles
h=ishold;
hold on;
s1=circle(center1+ cfg.position,r1,cfg.resolution,cfg.style{1});
s2=circle(center2+ cfg.position,r2,cfg.resolution,cfg.style{2});
s3=circle(center3+ cfg.position,r3,cfg.resolution,cfg.style{3});
%-- returns hold
if ~h, hold off;end 

%-- return answers
cfg.radius  = [r1 r2 r3];
cfg.centers = [center1 ;center2; center3]+ repmat(cfg.position,3,1);
cfg.handles = [s1 s2 s3];