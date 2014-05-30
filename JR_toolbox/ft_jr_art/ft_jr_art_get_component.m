function [Y cfg] = ft_jr_art_get_component(X,cfg)

if ~isfield(cfg,'rm_badchan'),  error('missing rm_badchan');    end
if ~isfield(cfg,'trl_sel'),     error('missing trl_sel');       end
if ~isfield(cfg,'void'),        cfg.void = true;                end
if ~isfield(cfg,'threshold'),   cfg.threshold = .15;            end

%% Define generic functions
%------ compatibility
% if ~exist('trimmean','file'),   trimmean = @(x,percent,flag,dim) nanmean(x,dim); end
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,15) ' ' x ' ' repmat('-',1,55-length(x))]);
else
    print = @(x) false;
end

%% average across trials
print('1. Mean artefact');
[n_chans n_time] = size(X{1});
n_trials = length(cfg.trl_sel);
y = reshape(cell2mat(X(cfg.trl_sel)),[n_chans n_time, n_trials]);
Y.avg = squeeze(trimmean(y,90,'round',3));

%% 2 remove bad chan
print(['2. Removed ' num2str(length(cfg.rm_badchan)) ' bad channel(s)']);
Y.avg(cfg.rm_badchan,:)  = 0;

%% 3. Compute artefact pca ---------------------------
print('3. Computes PCA across channels');
[Y.coeff...
    Y.score...
    Y.latent]= princomp(Y.avg');
% latent into explained variance
Y.latent = Y.latent ./ sum(Y.latent);

%% 4. Find relevant component
print('4. Find relevant components');
Y.rm_latent           = find(Y.latent >= cfg.threshold) ;
