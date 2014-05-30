function d = dprime(hit,fa,miss,cr)
% dprime(hit,fa,miss,cr)

% data.hit
% data.miss
% data.cr
% data.fa
if nargin == 1, 
    data = hit;
else
    data.hit = hit;
    data.fa = fa;
    data.miss = miss;
    data.cr = cr;
end
HR      = data.hit ./ (data.hit + data.miss);
FAR     = data.fa ./ (data.fa + data.cr);
d       = norminv(HR) - norminv(FAR) ;

return