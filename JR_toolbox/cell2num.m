function num = cell2num(cell)

for i = 1:numel(cell)
%     display(i);
    if isempty(cell2mat(cell(i)))
        num(i) = NaN;
    else
        num(i) = cell2mat(cell(i));
    end
end

num = reshape(num,size(cell));
return