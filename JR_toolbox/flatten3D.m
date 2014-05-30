function data = flatten3D(data,options)
% data = 3D vector: instance (lines) * x (column 1) * y (column 2) * z
% (column 3)
% return 2D vector
if nargin == 1
    options.tmp = [];
end
options.factor = getoption(options,'factor', .98);

%-- normalize
data(:,1) = (data(:,1) - mean(data(:,1)))./ std(data(:,1));
data(:,2) = (data(:,2) - mean(data(:,2)))./ std(data(:,2));
data(:,3) = (data(:,3) - mean(data(:,3)))./ std(data(:,3));
%-- plot with max on z
data(:,3) = data(:,3) - max(data(:,3));

%-- change 
data(:,1) = data(:,1) - data(:,1).* (data(:,3) .^ options.factor);
data(:,2) = data(:,2) - data(:,2).* (data(:,3) .^ options.factor);
data(:,3) = [];
end