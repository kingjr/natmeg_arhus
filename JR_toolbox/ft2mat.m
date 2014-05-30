function [data info] = ft2mat(data)
info = rmfield(data,'trial');
data = reshape(cell2mat(data.trial),size(data.trial{1},1), size(data.trial{1},2), []);
return