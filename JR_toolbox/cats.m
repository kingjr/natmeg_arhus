function m = cats(structure,field)
% m = cats(structure,field)
% extracts data from structure and reshape size properly 
if size(structure,1) == 1 && size(structure,2) > 1
    structure = structure';
end
eval(['m = [structure.' field '];'])
eval(['s = size(structure(1).' field ');']);
S = size(structure);
% m = reshape(m,[s S]);
m = reshape(m,[s(1:2) S s(3:end)]);
if length(size(m))>3
    m = permute(m,[1 2 (3+length(S)):length(size(m)) 3:(3+length(S)-1)]);
end