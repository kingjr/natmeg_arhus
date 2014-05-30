function binsave(filename, matrix)
% binsave(filename, matrix)
    fid = fopen( [filename '_'],'w');
    fwrite(fid, matrix, 'float');
    fclose(fid);
    movefile([filename '_'], filename);
end