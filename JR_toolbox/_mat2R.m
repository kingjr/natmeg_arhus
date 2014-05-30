function [vars log] = mat2R(X,command,cfg)
% [vars log] = mat2R(X,command [,cfg])
%--------------------------------------------------------------------------
% input:
%       X:              matrix or cell of n observation by m conditions
%       command:        R command to be executed
%       cfg.names:      names of each conditions
%       cfg.iv:         vector indidcating independent variables [2:end]
%       cfg.num:        vector indidcating numerical variables [1]
%       cfg.data_file:  file in which data will be saved
%       cfg.save_file:  file in which results will be saved
%       cfg.script_file:file in which script will be saved
%       cfg.log_file:   file in which R log will be save
%       cfg.delete_file:delete files at the end [true]
% ouput:
%       txt:            R txt output
%       log:            R log output
% 
% example 1:
%       cfg.names   = {'DV', 'IV1', 'IV2', 'rnd_factor'};
%       cfg.vars    = {'m'}; % retrieve m
%       result      = mat2R(cat(2,DV,IV1,IV2,rnd_factor),...
%           'm<-coef(summary(lmer(DV ~ IV1 * IV2 + (1|rnd_factor))))',cfg);
%
% example 2: chi square:
%       cfg         = [];
%       cfg.names   = {'A', 'B', 'C'};
%       cfg.num     = 1:3;
%       cfg.vars    = {'m'};
% vars = mat2R(data,'m<-chisq.test(rbind(A,B,C))',cfg);
% 
%--------------------------------------------------------------------------
% author: Jean-Remi King: jeanremi.king+matlab [at] gmail.com
%--------------------------------------------------------------------------
% 2012 11 14: add example, fix retrieve variable
% 2012 10 11: fix matrix 2 cell convertions
% 2012 09 26: main script
%--------------------------------------------------------------------------
%% Parameters
if ~exist('cfg','var'),         cfg             = [];                       end
if ~isfield(cfg,'names'),       cfg.names       = str2cell(1:size(X,2));    end
if ~isfield(cfg,'data_file'),   cfg.data_file   = [pwd '/_data.csv'];       end
if ~isfield(cfg,'script_file'), cfg.script_file = [pwd '/_script.R'];       end
if ~isfield(cfg,'save_file'),   cfg.save_file   = [pwd '/data.txt'];        end
if ~isfield(cfg,'save_mat'),    cfg.save_mat    = [pwd '/save.mat'];        end
if ~isfield(cfg,'log_file'),    cfg.log_file    = [pwd '/log.txt'];         end
if ~isfield(cfg,'iv'),          cfg.iv          = 2:size(X,2);              end
if ~isfield(cfg,'num'),         cfg.num         = 1;                        end
if ~isfield(cfg,'package'),     cfg.package     = {'arm', 'lme4', 'R2WinBUGS', 'coda', 'Rlab'};end
if ~isfield(cfg,'vars'),        cfg.vars        = {};                       end
if ~isfield(cfg,'delete_file'), cfg.delete_file = true;                     end

%% Format data
if isnumeric(X), X = num2cell(X); end

%% export data to CSV
%-- write header
fid=fopen(cfg.data_file,'wt');
for column = 1:size(cfg.names,2)
    if column == size(cfg.names,2)
        fprintf(fid,'%s',cfg.names{column});
    else
    fprintf(fid,'%s,',cfg.names{column});
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
fid=fopen(cfg.script_file,'wt');
for ii = 1:length(cfg.package)
    fprintf(fid,['library(' cfg.package{ii} ')\n']);
end
if ~isempty(cfg.vars)
    fprintf(fid,'library(R.matlab)\n');
end
fprintf(fid,['setwd(''' pwd ''');\n']);
fprintf(fid,['data=read.table(''' cfg.data_file ''', sep='','', header=T)\n']);
fprintf(fid,['for (i in c(' regexprep(num2str(cfg.iv),'  ', ',') '))\n']);
fprintf(fid,'  data[,i]=as.factor(data[,i])\n');
fprintf(fid,['for (i in c(' regexprep(num2str(cfg.num),'  ', ',') '))\n']);
fprintf(fid,'  data[,i]=as.numeric(as.vector(data[,i]))\n');
fprintf(fid,'attach(data)\n');
fprintf(fid,['capture.output(' command ', file = "' cfg.save_file '")\n']);
if ~isempty(cfg.vars)
    txt = ['writeMat("' cfg.save_mat '"'];
    for ii = 1:length(cfg.vars)
        txt = cat(2,txt, [', ' cfg.vars{ii} '=' cfg.vars{ii}]);
    end
    txt = cat(2,txt, ')');
    fprintf(fid,txt);
end
fclose(fid);

%% Prepare results
if exist(cfg.save_file,'file')
    delete(cfg.save_file);
end

%% Execute R
[status slog] = system(['R CMD BATCH ' cfg.script_file ' ' cfg.log_file]);

%% Retrieve results
log=char(textread(cfg.log_file,'%s','delimiter','\n'));
if ~exist(cfg.save_file,'file')
    txt=log;
    disp(slog);
    disp(log);
else
    txt=char(textread(cfg.save_file,'%s','delimiter','\n'));
end
if ~isempty(cfg.vars) 
    if exist(cfg.save_mat,'file')
        vars = load(cfg.save_mat);
    else
        disp(slog);
        disp(log);        
        vars = [];
    end
else
    vars = txt;
end
%% clean up
if cfg.delete_file
    delete(cfg.data_file);
    delete(cfg.script_file);
    delete(cfg.save_file);
    delete(cfg.log_file);
end