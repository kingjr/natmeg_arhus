function current = versioning_field(current,previous)
% current = versioning_field(current,previous)
cfields = fieldnames(current);
pfields = fieldnames(previous);
current.previous = previous;
for f = 1:length(pfields)
    if ~ismember(pfields{f},cfields)
        eval(['current.' pfields{f} '= previous.' pfields{f} ';']);
    end
end

end