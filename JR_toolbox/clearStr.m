function string = clearStr(string)
string = strrep(string, 'é', 'e');
string = strrep(string, 'ê', 'e');
string = strrep(string, 'è', 'e');
string = strrep(string, 'à', 'a');
string = strrep(string, 'â', 'a');
string = strrep(string, 'ç', 'c');
string = strrep(string, 'ù', 'u');
string = strrep(string, 'û', 'u');
string = strrep(string, 'ô', 'o');
string = strrep(string, 'î', 'i');
return