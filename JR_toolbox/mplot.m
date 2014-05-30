function mplot(curve)
if length(size(curve)) == 3
for ii = 1:size(curve,1)
    plot(curve(ii,:,1),curve(ii,:,2),'linewidth', 2.5, 'color', 'k');
end
elseif length(size(curve)) == 2
    plot(curve(:,1),curve(:,2), 'linewidth', 2.5, 'color', 'k');
end
end