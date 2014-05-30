function name = ttl2name(nb,data)

% returns the name value given a ttl:
% name = ttl2name(nb,data)
% data.names    => names from name2ttl
% data.dec      => ttl values from name2ttl


%-- build dimensions to find correct location
dims = '[';
for ii = 1:length(data.names)
    dims = [dims 'd' num2str(ii) ','];
end
dims(end) = ']';

%-- find value in names
eval([ dims '  = ind2sub(size(data.dec),find(data.dec==nb));']);

%-- return values for each dimension
name = '';
for dim = 1:length(data.names)
    name = [name ',' data.names{dim}{eval(['d' num2str(dim)])}];
end
name(1) = ''; % remove first coma
return