function d = cmpstr(str1,str2)
% compares string1 and string2, and return the differences
% str1 and str2 have to have identical length
index = find(str1~=str2);
d = [str1(index) str2(index)];