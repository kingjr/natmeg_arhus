function x = npol2ncart(r,th)
% npol2ncart Transforms n-dimensional polar coordinates to n-dimensional
% Cartesian coordinates.
%--------------------------------------------------------------------------
% x = npol2ncart(radius,theta)
% input
%       radius: vector of radius (each raw is a new point)
%       angles: matrix of angles (each raw is a point, each column a
%       dimension)
% output
%       x: matrix of cartesian coordinates (each raw is a point, each
%       column a dimension)
%--------------------------------------------------------------------------
% example:
%       radius = [1:1000]';
%       angles = [ones(1,1000); linspace(0,10*pi,1000); linspace(0,10*pi,1000)]';
%       x = npol2ncart(radius,angles);
%       scatter3(x(:,1),x(:,2),x(:,3));
%--------------------------------------------------------------------------
% see http://en.wikipedia.org/wiki/N-sphere for explanation
%--------------------------------------------------------------------------
% (c) Jean-Rémi King, jeanremi.king+matlab@gmail.com
%--------------------------------------------------------------------------
if size(th,1) ~= size(r,1), error('radius and thetas should have compatible dimensionalities'); end
n_points    = size(th,1);
n_dim       = size(th,2);
x           = NaN(n_points, n_dim); 
s           = sin(th);
c           = cos(th);
for dim = 1:(n_dim -1)
    x(:,dim) = r.* prod(s(:,1:(dim-1)),2).*c(:,dim);
end
x(:,n_dim)  = r.* prod(s(:,1:(n_dim-1)),2);