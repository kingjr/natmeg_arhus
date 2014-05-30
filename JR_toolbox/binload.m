function m = binload(filename, dims)
% m = binload(filename, dims)
    fid = fopen(filename,'r');
    a = fread(fid, 'float');
    fclose(fid);
    try
        m = reshape(a, dims);
    catch
        warning(['Could not reshape with these dimensions!']);
        m = a;
    end
end