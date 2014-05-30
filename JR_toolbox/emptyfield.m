
function out = emptyfield(struct, field, fill)
% out = emptyfield(struct, field, [fill])
% checks for empty field and optionally replace by a given var
if nargin == 2, fill = []; end
% convert to cell
eval(['X = {struct.' field '};']);
% find empty array
index = find(cell2mat(cellfun(@(x) isempty(x), X, 'uniformoutput', false)));
if nargin == 2, 
    out = index;
else
    for ii = 1:length(index)
        eval(['struct(index(ii)).' field '=fill;']);
    end
    out = struct;
end
    
