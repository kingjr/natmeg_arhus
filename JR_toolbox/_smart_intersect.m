function output = smart_intersect(cells)
output = cells{1};
for ii = 2:length(cells)
    output = intersect(output,cells{ii});
    
end

return