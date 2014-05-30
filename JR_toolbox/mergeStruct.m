function strct = mergeStruct(strct1,strct2)
% strct = mergeStruct(strct1,strct2)
% merges two structure into one 
% strct1(n), strct2(m) => strct(n+m)
strct = strct1;
for ii = 1:length(strct2)
    n = fieldnames(strct2(ii));
    for f = 1:length(n)
        eval(['strct(length(strct1)+ii).' n{f} '= strct2(ii).' n{f} ';']);
    end
end