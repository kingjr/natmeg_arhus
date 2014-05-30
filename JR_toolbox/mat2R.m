function out = mat2R(X,command,varargin)
% out = mat2R(X,command [,cfg])
%--------------------------------------------------------------------------
% input:
%       X:              matrix or cell of n observation by m conditions
%       command:        R command to be executed
%       names:      names of each conditions
%       iv:         vector indidcating independent variables [2:end]
%       num:        vector indidcating numerical variables [1]
%       data_file:  file in which data will be saved
%       save_file:  file in which results will be saved
%       script_file:file in which script will be saved
%       log_file:   file in which R log will be save
%       delete_file:delete files at the end [true]
% ouput:
%       txt:            R txt output
%       log:            R log output
% 
% example 1:
%       names   = {'DV', 'IV1', 'IV2', 'rnd_factor'};
%       vars    = {'m'}; % retrieve m
%       result      = mat2R(cat(2,DV,IV1,IV2,rnd_factor),...
%           'm<-coef(summary(lmer(DV ~ IV1 * IV2 + (1|rnd_factor))))',cfg);
%
% example 2: chi square:
%       cfg         = [];
%       names   = {'A', 'B', 'C'};
%       num     = 1:3;
%       vars    = {'m'};
% vars = mat2R(data,'m<-chisq.test(rbind(A,B,C))',cfg);
%
% example 3: anova
%       summary(aov(X ~ state * delay + Error(uid)))
% 
%--------------------------------------------------------------------------
% author: Jean-Remi King: jeanremi.king+matlab [at] gmail.com
%--------------------------------------------------------------------------
% 2012 11 14: add example, fix retrieve variable
% 2012 10 11: fix matrix 2 cell convertions
% 2012 09 26: main script
%--------------------------------------------------------------------------
if nargin == 2
    varargin = {};
end
for ii = 1:2:length(varargin)
    eval([varargin{ii} '=varargin{ii+1};']);
end
   
%% Parameters
if ~exist('names', 'var'),       names       = str2cell(1:size(X,2));    end
if ~exist('data_file', 'var'),   data_file   = [pwd '/_data.csv'];       end
if ~exist('script_file', 'var'), script_file = [pwd '/_script.R'];       end
if ~exist('save_file', 'var'),   save_file   = [pwd '/data.txt'];        end
if ~exist('save_mat', 'var'),    save_mat    = [pwd '/save.mat'];        end
if ~exist('log_file', 'var'),    log_file    = [pwd '/log.txt'];         end
if ~exist('iv', 'var'),          iv          = 2:size(X,2);              end
if ~exist('num', 'var'),         num         = 1;                        end
if ~exist('package', 'var'),     package     = {'arm', 'lme4', 'R2WinBUGS', 'coda', 'Rlab'};end
if ~exist('vars', 'var'),        vars        = {};                       end
if ~exist('delete_file', 'var'), delete_file = true;                     end

%% Format data
if isnumeric(X), X = num2cell(X); end

%% export data to CSV
%-- write header
fid=fopen(data_file,'wt');
for column = 1:size(names,2)
    if column == size(names,2)
        fprintf(fid,'%s',names{column});
    else
    fprintf(fid,'%s,',names{column});
    end
end
%-- write data
fprintf(fid,'\n');
for raw=1:size(X,1)
    for column = 1:size(X,2)
        if isnumeric(X{raw,column}) || islogical(X{raw,column}), X{raw,column}=num2str(X{raw,column}); end
        if column == size(X,2)
            fprintf(fid,'%s',X{raw,column});
        else
            fprintf(fid,'%s,',X{raw,column});
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

%% Create R file
fid=fopen(script_file,'wt');
for ii = 1:length(package)
    fprintf(fid,['library(' package{ii} ')\n']);
end
if ~isempty(vars)
    fprintf(fid,'library(R.matlab)\n');
end
fprintf(fid,['setwd(''' pwd ''');\n']);
fprintf(fid,['data=read.table(''' data_file ''', sep='','', header=T)\n']);
fprintf(fid,['for (i in c(' regexprep(num2str(iv),'  ', ',') '))\n']);
fprintf(fid,'  data[,i]=as.factor(data[,i])\n');
fprintf(fid,['for (i in c(' regexprep(num2str(num),'  ', ',') '))\n']);
fprintf(fid,'  data[,i]=as.numeric(as.vector(data[,i]))\n');
fprintf(fid,'attach(data)\n');
fprintf(fid,['capture.output(' command ', file = "' save_file '")\n']);
if ~isempty(vars)
    txt = ['writeMat("' save_mat '"'];
    for ii = 1:length(vars)
        txt = cat(2,txt, [', ' vars{ii} '=' vars{ii}]);
    end
    txt = cat(2,txt, ')');
    fprintf(fid,txt);
end
fclose(fid);

%% Prepare results
if exist(save_file,'file')
    delete(save_file);
end

%% Execute R
[out.status out.slog] = system(['R CMD BATCH ' script_file ' ' log_file]);

%% Retrieve results
out.log=char(textread(log_file,'%s','delimiter','\n'));
if exist(save_file,'file')
    out.txt=char(textread(save_file,'%s','delimiter','\n'));
end
if ~isempty(vars) 
    if exist(save_mat,'file')
        out.vars = load(save_mat);
    end
end
%% clean up
if delete_file
    delete(data_file);
    delete(script_file);
    delete(save_file);
    delete(log_file);
end