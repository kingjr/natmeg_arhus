function csvwritecell(filename,data,varargin)
% csvwritecell(filename,data,varargin)
% Jean-Rémi King
if nargin == 2, varargin = {}; end
for ii = 1:2:length(varargin), eval([varargin{ii} '=varargin{ii+1};']); end
fid=fopen(filename,'wt');
[nl nc] = size(data);
for raw=1:nl
    for column = 1:nc
        if isnumeric(data{raw,column}) || islogical(data{raw,column}), 
            data{raw,column}=num2str(data{raw,column}); 
        end
        if column == nc
            fprintf(fid,'%s\n',data{raw,column});
        else
            fprintf(fid,'%s,',data{raw,column});
        end
    end
end
fclose(fid);
