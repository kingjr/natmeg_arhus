function s = sz(x,dims)
% s = sz(x[,dims])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of matrix and gives multiple specific dimensions
% e.g. 
%       sz(rand(2,3,4),[3 4])
%       sz(rand(2,3,4),'3:end');
%       sz(rand(2,3,4),[1 4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jeanremi [dot] king [at] gmail [dot] com

if nargin == 1
    s = size(x);
else
    s = size(x);
    if ischar(dims)
        eval(['s = s(' dims ');']);
    else
        s = s(dims);
    end
end