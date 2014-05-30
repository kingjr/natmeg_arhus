function output = summary(data)

f= '';
var = fieldnames(data);

for field = 1:length(var)
    %-- get unique values in fields
    cdata = struct2cell(data);
    cdata = squeeze(cdata(strcmp(fieldnames(data), field),:,:));

    f = ['for var{field} = 1: f end']
end


return