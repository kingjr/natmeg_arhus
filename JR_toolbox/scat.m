function output = scat(dim,input)
% scat(dim,{input})
output = [];
for ii = 1:length(input)
    output = cat(dim,output,input{ii});
end
return