function dec_codes = name2ttl(data)
% transform a list of cells into ttl codes
% dec_codes = name2ttl(data)
% eg. data: {{'conditionA1', 'conditionA2},
% {'conditionB1','conditionB2','conditionB3'}}

[str1 str2] = deal('');
for cond = 1:length(data)
    str1 = [str1 ',' num2str(length(data{cond}))];
    str2 = [str2 '*' num2str(length(data{cond}))];
end
eval(['dec_codes = reshape(1:(' str2(2:end) ')' str1 ');']);
if dec_codes(end) > 31, warning('more than 32 values!');end
return