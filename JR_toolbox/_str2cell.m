function cells = str2cell(txt)
for ii = 1:length(txt)
    cells{ii} = txt(ii);
end
return