function all_nb = embed_pie_sum(data,level,all_nb,current_level)

if nargin == 2
    all_nb = [];
    current_level = 1;
end

    
if iscell(data(1))  
    for index = 1:length(data)
        all_nb = embed_pie_sum(data{index},level,all_nb,current_level+1);
    end
else
    if current_level >= level
        all_nb = cat(2,all_nb, data);
    end
end
return