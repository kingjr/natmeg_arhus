function out = ft_neighbor_nD(pos,labels,distance_threshold)
% find neighbouring channels in a N dimensional space

%% compute distances 
nchan = size(pos,1);
distances = NaN(nchan,nchan);
for c1 = 1:(nchan-1)
    for c2 = (c1+1):nchan
        distances(c1,c2) = sqrt(sum((pos(c1,:)-pos(c2,:)).^2));
        distances(c2,c1) = distances(c1,c2);
    end
end

%% find neighbor
for c = nchan:-1:1
    neighbors           = distances(c1,:)<distance_threshold;
    out(c).label        = labels{c};
    out(c).neighblabel  = labels(neighbors);
    n(c)                = sum(neighbors);
end

%% output 
disp(['Using the ' num2str(size(pos,2)) '-D layout to determine the neighbours']);
disp(['There are on average ' num2str(nanmean(n)) ' neighbours per channel']);