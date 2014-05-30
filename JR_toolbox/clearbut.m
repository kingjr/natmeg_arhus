vars = who;
ii = 1;
keepVars{length(keepVars)+1} = 'ii';
keepVars{length(keepVars)+1} = 'keepVars';
while length(who) > length(keepVars)
    vars = who;
    if ~ismember(vars{ii}, keepVars)
        eval(['clear ' vars{ii} ';']);
    else
        ii = ii + 1;
    end
end
clear ii keepVars
