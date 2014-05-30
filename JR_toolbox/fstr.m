function index = fstr(db,s) 
% fast multiple string search in cell array
% JeanRémi King
if ischar(s), s = {s}; end
index = cell2mat(cellfun(@(x) find(ismember(db,x)), s, 'uniformoutput', false));