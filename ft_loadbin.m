function data = ft_loadbin(file_ft,file_dat)
load(file_ft,'data');
dat = binload(file_dat,data.dims);
data.trial = {};
for trial = data.dims(3):-1:1
    data.trial{trial} = dat(:,:,trial);
    dat(:,:,trial) = 0;
end
end

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