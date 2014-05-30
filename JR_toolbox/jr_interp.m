function curve=jr_interp(prop,curve1,curve2,sizes)
if nargin == 3
    sizes = [1 1;1 1];
end
%-- remove NaN 
curve1 = curve1(intersect(find(~isnan(curve1(:,1))),find(~isnan(curve1(:,2)))),:);
curve2 = curve2(intersect(find(~isnan(curve2(:,1))),find(~isnan(curve2(:,2)))),:);

%-- remove duplicates
dupli1 = NaN;dupli2 = NaN;
while ~isempty([dupli1;dupli2])
    dupli1 = intersect(find(curve1(1:end-1,1) == curve1(2:end,1)),find(curve1(1:end-1,2) == curve1(2:end,2)));
    dupli2 = intersect(find(curve2(1:end-1,1) == curve2(2:end,1)),find(curve2(1:end-1,2) == curve2(2:end,2)));
    curve1(dupli1(1:2:end),:) = [];
    curve2(dupli2(1:2:end),:) = [];
end

%-- check that curve have the same number of dots;
if length(curve1) ~= length(curve2)
    curve1 = interparc(max([length(curve1),length(curve2)]),curve1(:,1),curve1(:,2));
    curve2 = interparc(max([length(curve1),length(curve2)]),curve2(:,1),curve2(:,2));
end

%-- interpolation
curve = [];
for s = prop
    %-- linear interpolation
    curve(end+1,:,1:2) = [curve1(:,1).*(1-s) + curve2(:,1).*s,curve1(:,2).*(1-s)+curve2(:,2).*s];
end
curve = squeeze(curve);