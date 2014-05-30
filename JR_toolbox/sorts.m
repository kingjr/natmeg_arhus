function new_struct = sorts(old_struct, value,option)
if nargin == 2
    option = 'ascend';
end
[unused, order] = eval(['sort([old_struct(:).' value '], ''' option ''')']);
new_struct = old_struct(order); 