function cfg = getscripts(cfg)
% cfg = getscripts([cfg])
% read and output all files in a variable to keep track of file changes
% by default:
%       cfg.path = pwd
%       cfg.filetype = '*.m'
% JeanRémi King, jeanremi.king@gmail.com

if nargin == 0,
    cfg = [];
end
if ~isfield(cfg, 'path'), cfg.path = [pwd '/']; end
if ~isfield(cfg, 'filetype'), cfg.filetype = '*.m'; end

%-- read all script in folder
files = dir([cfg.path cfg.filetype ]);
for f = 1:length(files) 
    fid = fopen([cfg.path files(f).name], 'r');
    try
        %eval(['cfg.' files(f).name(1:(end-(length(cfg.filetype)-1))) ' = fread(fid,''uint8=>char'')'';']);
        eval(['cfg.' files(f).name ' = fread(fid,''uint8=>char'')'';']);
    catch
        % deal with bad naming
        eval(['cfg.file_' num2str(f) ' = fread(fid,''uint8=>char'')'';']);
    end
    fclose(fid);
end
