function addpathr(path)
% addpathr(path)
% add path with all subfolders
% works on unix systems 
% (c) JeanRemi King 2011

if strcmp(path(end),'/'), path(end) = [];end    % to avoid bugs, remove last '/'
addpath(path);
files = dir(path);                              % list containing elements
for file = 3:length(files)                      % avoid ../ and ./ folders
    if files(file).isdir
        addpathr([path '/' files(file).name]);  % add folders
    end
end
