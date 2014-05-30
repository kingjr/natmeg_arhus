function index = cfind(string,carray)
% index = cfind(string,carray)
% finds strings in a cell array
% (c) JeanRemi King, 2012.
% jeanremi.king [at] gmail [dot] com
scfind =@(string,carray) find(cell2mat(cellfun(@(x) ~isempty(strfind(x,string)),carray,'UniformOutput',false))); 
if iscell(string)
    index = cell(size(string));
    for c = 1:numel(string)
        index{c} = scfind(string{c},carray);
    end
else
    index = scfind(string,carray);
end
