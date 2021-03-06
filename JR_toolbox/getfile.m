function out = getfile(files)
% cfg = getscripts({filenames})
% read and output all files in a variable to keep track of file changes
% JeanRémi King, jeanremi.king@gmail.com
out = [];
if ~iscell(files), 
    if strfind(files, '*') || strcmp(files(end), '/')
        s = findstr(files, '/');
        path = files(1:s(end));
        files = dir(files);
        files = {files.name};
        files = cellfun(@(x) [path x], files, 'UniformOutput', false);
    else
    files = {files}; 
    end
end

%-- read all script in folder
for f = 1:length(files) 
    fid = fopen(files{f}, 'r');
    s = findstr(files{f}, '/');
    details = dir(files{f});
    try
        out(f).file = fread(fid,'uint8=>char')';
        fclose(fid);
    catch
        out(f).file = 'ERROR INVALID ENCODING';
    end
        
    out(f).filename = details.name;
    out(f).date = details.date;
    out(f).datenum = details.datenum;
    out(f).bytes = details.bytes;
    out(f).path = files{f}(1:s(end));
    
end
