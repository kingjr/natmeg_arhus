function  ft_jr_art_plot(component,cfg)
% ft_jr_art_plot(component[,cfg])
%
% plot artefact detection and correction results
%
% requires: ft_jr_art.m
%--------------------------------------------------------------------------
% Jean-Rémi King
%--------------------------------------------------------------------------


if nargin == 1
    cfg = [];
end

%
if ~isfield(cfg,'plot_art_detect'),                     cfg.plot_art_detect = true;     end
if ~isfield(cfg,'plot_mean_art'),                       cfg.plot_mean_art= true;        end
if ~isfield(cfg,'plot_mean_art_topo'),                  cfg.plot_mean_art_topo= true;   end
if ~isfield(cfg,'plot_latent_corr'),                    cfg.plot_latent_corr= true;     end
if ~isfield(cfg,'plot_mean_art_topo_gen'),              cfg.plot_mean_art_topo_gen= true;   end
if ~isfield(cfg,'plot_mean_art_topo_gen_corrected'),    cfg.plot_mean_art_topo_gen_corrected= true;   end
if ~isfield(cfg,'plot_selected_topo'),                  cfg.plot_selected_topo= true;   end

n_subplots = cfg.plot_art_detect +....
    cfg.plot_mean_art + ...
    cfg.plot_mean_art_topo + ...
    cfg.plot_latent_corr + ...
    cfg.plot_mean_art_topo_gen + ...
    cfg.plot_mean_art_topo_gen_corrected+...
    cfg.plot_selected_topo;
sb = 0;


% plot layout artefact
if ~isfield(component.cfg, 'layout'), component.cfg.layout = cell(length(component.cfg.chantypes),1); end

% automatically determines scale
if ~isfield(cfg,'scales'),
    if length(component.cfg.chantypes) == 1 && strcmpi(component.label{1}(1:3), 'EEG')
        cfg.scale = [-1 1] .* 10^-3;
    elseif length(component.cfg.chantypes) == 1 && strcmpi(component.label{1}(1:3), 'MEG')
        cfg.scale = [-1 1] .* 10^-11;
    elseif length(component.cfg.chantypes) == 2 && strcmpi(component.label{1}(1:3), 'MEG')
        cfg.scale = [[-5 5] .* 10^-12; [-1 1] .* 10^-12];
    end
end

clf;set(gcf,'name',component.cfg.dataset,'color','w');
%-- artefact finding
if cfg.plot_art_detect
    sb =sb+1;
    subplot(n_subplots,1,sb);cla;hold on;
    time = (1:size(component.trials.artchan,2))/component.trials.cfg.fsample;
    plot(time,component.trials.artchan,'b');
    scatter(time(component.trials.all_trl(component.trials.trl_sel.train)), ...
        component.trials.cfg.sample_threshold*ones(length(component.trials.trl_sel.train),1),'*r');
    scatter(time(component.trials.all_trl(component.trials.trl_sel.test)),...
        component.trials.cfg.sample_threshold*ones(length(component.trials.trl_sel.test),1),'*g');
  
    axis([min(time) max(time) ylim]);
    title('artefact channel and artefact detection');box off;
end
for chantype = 1:length(component.cfg.chantypes)
    sb = cfg.plot_art_detect-1;
    %-- average artefact
    if cfg.plot_mean_art
        sb = sb+1;
        subplot(n_subplots,length(component.cfg.chantypes),length(component.cfg.chantypes)*sb+chantype);cla;
        plot(component.trials.time, component.trials.artchan_mean,'b');
        axis([min(component.trials.time) max(component.trials.time) ylim]);
        title(['mean artefeact component']);box off;
    end
    %-- average ERP
    if cfg.plot_mean_art_topo
        sb=sb+1;
        subplot(n_subplots,length(component.cfg.chantypes),length(component.cfg.chantypes)*sb+chantype);cla;
        imagesc(component.trials.time,component.cfg.chantypes{chantype},component.comp.avg(component.cfg.chantypes{chantype},:),cfg.scale(chantype,:));%normal
        title(['ERP ' num2str(chantype)]);box off;
    end
    %-- correlation & latent values
    if cfg.plot_latent_corr
        sb = sb+1;
        subplot(n_subplots,length(component.cfg.chantypes),length(component.cfg.chantypes)*sb+chantype);cla; hold on;
        scatter(abs(component.corr.R{chantype}),component.comp.pca(chantype).latent, '+b');
        scatter(abs(component.corr.R{chantype}(1:size(component.correct.remove{chantype},1))),component.comp.pca(chantype).latent(1:size(component.correct.remove{chantype},1)), 'r', 'filled');
        axis([0 1 0 1]);box off;xlabel('R'); ylabel('Latent');
    end
    
    %-- original independent
    if cfg.plot_mean_art_topo_gen
        sb = sb+1;
        subplot(n_subplots,length(component.cfg.chantypes),length(component.cfg.chantypes)*sb+chantype);cla;
        imagesc(component.trials.time,...
            component.cfg.chantypes{chantype},....
            component.correct.avg_dirty{chantype},cfg.scale(chantype,:));
        title(['Corrected ERP ' num2str(chantype)]);box off;
    end
    %-- correction
    if cfg.plot_mean_art_topo_gen_corrected
        sb = sb+1;
        subplot(n_subplots,length(component.cfg.chantypes),length(component.cfg.chantypes)*sb+chantype);cla;
         imagesc(component.trials.time,...
            component.cfg.chantypes{chantype},....
            component.correct.avg_clean{chantype},cfg.scale(chantype,:));
        title(['Corrected ERP ' num2str(chantype)]);box off;
    end
    %-- plot topo
    if cfg.plot_selected_topo
        sb = sb+1;
        if ~isempty(component.cfg.layout{chantype});
            % for each removed component
            for topo = 1:length(component.correct.remove{chantype})
                X = [component.comp.pca(chantype).coeff(:,component.correct.remove{chantype}(topo)); 0;0];
                cfg_plot = [];
                cfg_plot.layout = component.cfg.layout{chantype};
                cfg_plot.label = cfg_plot.layout.label;
                cfg_plot.zlim = [-3 3].*mad(X) + median(X);
                cfg_plot.marker= 'off';
                subplot(n_subplots,...
                    length(component.cfg.chantypes)*length(component.correct.remove{chantype}),...
                    length(component.cfg.chantypes)*length(component.correct.remove{chantype})*sb+...
                    (chantype-1)*length(component.correct.remove{chantype})+topo);cla;
                my_plot_topo(X,cfg_plot);
            end
        end
    end
end
return