function [Y cfg] = ft_jr_art_corr(comp,art,cfg)
% Compute components correlation with artefacted channel
if ~isfield(cfg,'void'),        cfg.void        = true;     end % display feedback
if ~isfield(cfg,'apply'),       cfg.apply       = true;     end % optional computation
if ~isfield(cfg,'threshold'),   cfg.threshold   = .05;      end
if ~isfield(cfg,'maxsample'),   cfg.maxsample   = 10000;    end % maximum number of sample to consider


%% Define generic functions
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,15) ' ' x ' ' repmat('-',1,55-length(x))]);
else
    print = @(x) false;
end

if cfg.apply
    print('1. Compute correlations');
    n_samples   = (size(comp,2)*size(comp,3));
    cfg.samples = 1:n_samples;
    cfg.samples = randperm(n_samples);
    % minimize number of samples used if necessary
    n_samples   = min(cfg.maxsample,n_samples);
    cfg.samples = cfg.samples(1:n_samples);
    %-- concatenate trials and samples
    comp        = reshape(comp,size(comp,1),[])';
    art         = reshape(art,1,[])';
    %-- compute correlation
    [Y.R Y.p]   = corr(comp(cfg.samples,:),art(cfg.samples));
else
    [Y.R Y.p] = deal(NaN(size(X,1),1));
end

%-- Find to-be-removed components
print('2. Find components correlated to artchan');
Y.threshold = cfg.threshold;
Y.rm_corr = find(Y.p <= cfg.threshold);
    