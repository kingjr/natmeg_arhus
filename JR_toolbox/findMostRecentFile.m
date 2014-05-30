function lf = findMostRecentFile(directory)

d = dir([directory '*.mat']);
[dx,dx]=sort([d.datenum],'descend');
lf=d(dx(1)).name;
return