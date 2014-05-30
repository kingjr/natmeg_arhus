function all=build_composite_img(colors,img,n_column,space)

img = imread(img);
all = [];
for jj = 1:ceil(length(colors)/n_column)
    index = ((jj-1)*n_column+1):((jj-1)*n_column+n_column);
    index(index>length(colors)) = [];
    line = newline(colors(index,:),img,space(1));
    if size(all,1)>0
        white = 255*ones(space(2),size(all,2),size(all,3));
        if size(line,2) ~= size(all,2)
            all = cat(1,all,white,cat(2,line, 255*ones(size(line,1),size(all,2)-size(line,2),size(all,3))));
        else
            all = cat(1,all,white,line);
        end
    else
        white = 255*ones(space(2),size(line,2),size(line,3));
        all = line;
    end
    
end



return

function line = newline(colors,img,space)
line = [];
white = 255*ones(size(img,1),space,size(img,3));
for ii = 1:length(colors)
    line = cat(2,line,white, newimg(colors(ii,:),img));
end

return

function img = newimg(color,img)
index = find(img<255);
c = repmat(permute(color,[1 3 2]),[size(img,1),size(img,2)]);
img(index) = double(255-img(index)).*double(c(index));
return