function string = ckName(string)
while ismember(' ',string)
    [t loc] = ismember(' ', string);
    string(loc) = '_';
end

