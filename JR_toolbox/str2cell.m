function cells = str2cell(txt)
for ii = 1:length(txt)
    if isnumeric(txt(ii))
        cells{ii} = num2str(txt(ii));
    else
        cells{ii} = txt(ii);
    end
end
return