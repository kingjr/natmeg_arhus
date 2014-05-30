function component = ft_jr_art_norm(component,cfg)
% component = ft_jr_art_norm(component,cfg)
% recompute PCA using options
if ~isfield(cfg,'norm'),      cfg.norm      = false; end % force PCA on a normalized across channel artefact
if ~isfield(cfg,'chantypes'), cfg.chantypes = component.cfg.chantypes; end
if ~isfield(cfg,'avg'),       cfg.avg       = component.avg; end
if ~isfield(cfg,'latent_thresh'),cfg.latent_thresh=component.cfg.latent_thresh;end % minimum variance explained by component (in %)

% 1. normalization
switch cfg.norm
    case 'zscore'
        disp('normalize zscore');
        cfg.avg = (cfg.avg - repmat(nanmean(cfg.avg,2),[1 size(cfg.avg,2)]))./repmat(nanstd(cfg.avg')',[1 size(cfg.avg,2)]);
    case 'medmad'
        disp('normalize medmad');
        cfg.avg = (cfg.avg - repmat(nanmedian(cfg.avg,2),[1 size(cfg.avg,2)]))./repmat(mad(cfg.avg')',[1 size(cfg.avg,2)]);
end
% 2. Recompute artefact pca
for chantype = 1:length(cfg.chantypes) % independently for each sensor
    disp(['channel type ' num2str(chantype)]);
    disp('computes pca ...');
    [component.pca(chantype).coeff...
        component.pca(chantype).score...
        component.pca(chantype).latent]= princomp(cfg.avg(cfg.chantypes{chantype},:)');
    % latent into explained variance
    component.pca(chantype).latent = component.pca(chantype).latent ./ sum(component.pca(chantype).latent);
    
    % reset correlation values
    [component.corr(chantype).R(:) component.corr(chantype).p(:)] = deal(0); 
    %-- Find to-be-removed components 
    component.rm_pca{chantype}           = find(component.pca(chantype).latent >= cfg.latent_thresh) ;
    component.rm_components{chantype}    = component.rm_pca{chantype};
    disp(['removed ' num2str(length(component.rm_components{chantype})) ' component(s)']);
    %-- Find to-be-kept components
    component.keep_components{chantype}  = setdiff(1:size(component.pca(chantype).coeff,1),component.rm_components{chantype});
    %-- Compute difference
    component.clear_comp{chantype}       = zeros(length(component.rm_components{chantype}),size(component.pca(chantype).coeff,1));
    component.clear_comp{chantype}       = cat(2,component.clear_comp{chantype}',component.pca(chantype).coeff(:,component.keep_components{chantype}))';
    
    %-- Compute correction on a different average artefact
    component.clear_art(cfg.chantypes{chantype},:)     = ...
        component.pca(chantype).coeff*...
        component.clear_comp{chantype}*...
        component.avg2(cfg.chantypes{chantype},:);
end
