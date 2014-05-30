function output = smart_intersect(cells)
% smart_intersect(struct)
% struct: {data, 'field', value)
output = cells{1};
for ii = 2:length(cells)
    output = intersect(output,cells{ii});
    
end

return