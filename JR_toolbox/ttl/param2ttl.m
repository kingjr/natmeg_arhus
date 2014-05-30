function dec = param2ttl(cfg,ttl)
% identify ttl digit from its structure
% dec = param2ttl(cfg,ttl)
fields = fieldnames(ttl.fields);
vector = []; % index
% identify, for each condition, the correct value
for f = 1:length(fields)
    if isnumeric(eval(['ttl.fields.' fields{f}]))
        vector(end+1) = find(eval(['ttl.fields.' fields{f} '== cfg.' fields{f}]));
    else
        vector(end+1) = find(ismember(...
            eval(['ttl.fields.' fields{f} ]),...
            eval(['cfg.' fields{f} ])));
    end
end
%-- report consequent ttl value from array
str = '';
for v = 1:length(vector)
    str = [str num2str(vector(v)) ',' ];
end
str(end) = []; % remove coma 
eval(['dec = ttl.dec(' str ');']);
return