function names = param2ttlnames(cfg)
% transform a structure into a series of names for ttl uses
% e.g. cfg.conditions1 = {'ala', 'bouni'}; cfg.condXX = 1:10;
fields = fieldnames(cfg);
names = '';
for field = 1:length(fields)
    names(end+1) = '{';
    for value = 1:length(eval(['cfg.' fields{field}]))
        if isnumeric(eval(['cfg.' fields{field}]))
            names = [names '''' fields{field} '_' num2str(eval(['cfg.' fields{field} '(value)'])) ''''];
        else
            names = [names '''' fields{field} '_' eval(['cfg.' fields{field} '{value}']) ''''];
        end
        names(end+1) = ',';
    end
    names(end:end+1) = '};';
end
names = eval(['{' names '};']);
return