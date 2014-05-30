function colorMatrix = colorCode(allLabels)
labels = unique(allLabels);
labelNb = length(labels); % number of different label
%colors = ([mod(3.*(0:labelNb-1)/labelNb, 1)' mod(3.*(0:labelNb-1)/labelNb, 2)' mod(3.*(0:labelNb-1)/labelNb, 3)']) .* 255/3;
colors = ...
    [0 0 3;0 3 0; 3 0 0;...
    0 1 2; 1 0 2; 0 2 1; 1 2 0; 2 0 1; 2 1 0; ...
    1 1 1; ...
    0 .5 2.5; .5 0 2.5; 0 2.5 .5; .5 2.5 0; 2.5 0 .5; 2.5 .5 0;] .* 255/3;
while size(colors,1) < labelNb
    colors(end+1,:) = [rand rand rand] .* 255;
end
    
for label = 1:labelNb % for each label
    index =findc(allLabels,labels(label));
    colorMatrix(index,1:3) = repmat(colors(label,:), length(index),1);
end
end