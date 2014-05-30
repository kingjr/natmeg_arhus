function output = smartCell2mat(data)
dim = length(size(data{1}));
for ii = 1:length(data)
    eval(['output(ii' repmat(',:',1,dim) ') = cell2mat(data(ii));'])
end

return