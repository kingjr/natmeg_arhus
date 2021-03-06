function [Y cfg] = ft_jr_art_generalize(X,cfg)

if ~isfield(cfg,'void'),        cfg.void     = true; end % display feedback
if ~isfield(cfg,'trl_sel'),     error('specify cfg.trl_sel'); end % display feedback
%% Define generic functions
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,15) ' ' x ' ' repmat('-',1,55-length(x))]);
else
    print = @(x) false;
end
%------ compatibility
% if ~exist('nanmean','file'),    nanmean= @(x,dim) mean(x,dim); warning('Could not find nanmean'); end
% if ~exist('trimmean','file'),   trimmean = @(x,percent,flag,dim) nanmean(x,dim); end


%% main
print('1. Average');
[n_chans n_time]    = size(X{1});
n_trials            = length(cfg.trl_sel);
Y.gen              = reshape(...
    cell2mat(X(cfg.trl_sel)),...
    [n_chans n_time, n_trials]);
Y.avg               = squeeze(trimmean(Y.gen,90,'round',3));

%-- rotate data into PC space
print('2. Rotate to component space');
Y.component         = reshape(cfg.coeff'*reshape(Y.gen,n_chans,[]),[size(cfg.coeff,1), n_time, n_trials]);
