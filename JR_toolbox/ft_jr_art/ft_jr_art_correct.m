function [Y cfg]= ft_jr_art_correct(X,cfg)

if ~isfield(cfg,'latent_thresh'),cfg.latent_thresh=.15;                     end % minimum variance explained by component (in %)
if ~isfield(cfg,'max_comp'),    cfg.max_comp    = 2;                        end % maximum principal components to be removed
if ~isfield(cfg,'void'),        cfg.void        = true;                     end % void
if ~isfield(cfg,'rm_corr'),     error('need cfg.rm_corr');                  end % removed components based on correlations R.
if ~isfield(cfg,'rm_latent'),   error('need cfg.rm_latent');                end % removed components based on latents

%% Define generic functions
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,15) ' ' x ' ' repmat('-',1,55-length(x))]);
else
    print = @(x) false;
end

%% Main
%-- Find to-be-removed components
Y.remove    = intersect(cfg.rm_corr, cfg.rm_latent);
print(['1. Removed ' num2str(length(Y.remove)) ' component(s)']);
%-- Find to-be-kept components
Y.keep  = setdiff(1:size(cfg.coeff,1),Y.remove);
%-- Compute difference
Y.clear       = zeros(length(Y.remove),size(cfg.coeff,1));
Y.clear       = cat(2,Y.clear',cfg.coeff(:,Y.keep))';

%-- Compute correction on average artefact
print('2. Correct on generalization data');
Y.avg_dirty = X;
Y.avg_clean = cfg.coeff*Y.clear*X;
return