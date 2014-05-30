function str = win2lin(str)

while ismember('\',str)
    [t pos] = ismember('\',str);
    str(pos) = '/';
end

return