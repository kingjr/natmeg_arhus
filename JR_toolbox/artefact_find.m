function output = artefact_find(cfg)
% output = ft_jr_ecgeog(cfg)
% -------------------------------------------------------------------------
% function that automatically removes the ECG artefact using a PCA on the
% average artefact.
% -------------------------------------------------------------------------
% inputs
%   - cfg.dataset           => path of dataset(s)       (can handle
%           multiple files in case you want to have a maximum number of art)
%   - cfg.plot              => 'yes' or 'no'            (default = 'yes')
%   - cfg.chantypes         => cells of chantypes       (default => all)
%   - cfg.ecgchan           => ecg channel              (default => read from header)
%   - cfg.dataformat        =>                          (default => 'neuromag')
%   - cfg.headerformat      =>                          (default => 'neuromag')
%   - cfg.prestim           => cut before ecg           (default => .2)
%   - cfg.poststim          => cut after ecg            (default => .5)
%   - cfg.dividetrial       => divide computation       (default => by 1)
%   - cfg.corr_thresh       => corr(PC,ecg) threshold   (default => .1)
% output
%   - output.chantypes      => same as above
%   - output.all_comp       => nchan x nchan matrix of artefact component
%   - output.clear_comp     => nchan x nchan matrix of non-artefacted comp.
%   - output.trial_nb       => number of artefact in each dataset file
%   - output.artefact_index => timing (in sample) of artefacts for each file
% -------------------------------------------------------------------------
% Requires fieldtrip 2011
% -------------------------------------------------------------------------
% (c) 2011 Jean-Rémi KING, all rights reserved
% jeanremi.king+matlab@gmail.com
% -------------------------------------------------------------------------

if ~isfield(cfg,'dataset'),     error('needs cfg.dataset');                 end
if ischar(cfg.dataset), cfg.dataset{1} = cfg.dataset;                       end
if ~isfield(cfg,'plot'),        cfg.plot        = 'yes';                    end
%-- read header
hdr                             = ft_read_header(cfg.dataset{1});
disp(hdr);
if ~isfield(cfg,'chantypes'),   cfg.chantypes   = {1:hdr.nChans};           end % cells dividing types of sensors (gradiometers, etc)
%-- finds ecg channels
if ~isfield(cfg,'artchan'),     cfg.artchan     = find(cell2mat(cellfun(@(x) ~isempty(findstr('EOG',x)), hdr.label, 'UniformOutput', false))); end
if isempty(cfg.artchan),        error('couldnt find cfg.artchan');          end
if ~isfield(cfg,'art_prestim'), cfg.art_prestim = .500;                     end
if ~isfield(cfg,'art_poststim'),cfg.art_poststim= 1.000;                    end
if ~isfield(cfg,'dividetrial'), cfg.dividetrial = 1;                        end % divide the computation by n trials for memory issue
if ~isfield(cfg,'corr_thresh'), cfg.corr_thresh = 3;                        end % correlation between ECG and PC, in STD
if ~isfield(cfg,'artnb'),       cfg.artnb       = 200;                      end % max number of artefacts to keep
if ~isfield(cfg,'threshold'),   cfg.threshold   = 3.5;                      end % in STD
if ~isfield(cfg,'max_comp'),    cfg.max_comp    = 3;                        end % maximum PC
cfg.trialfun                    = 'artefact_trialfun_thresh';               % get artefact beat as triggers

%-- find artefacts on each file and create trial structure
data_art                        = {};
data_raw                        = {};
artefact_index                  = {};
parfor file = 1:length(cfg.dataset)
    disp(['read file ' cfg.dataset{file} ' ...']);
    cfg_def                     = cfg;
    cfg_def.dataset             = cfg.dataset{file};
    cfg_def.trial               = 1:cfg.artnb;
    cfg_def.prestim             = cfg.art_prestim;
    cfg_def.poststim            = cfg.art_poststim;
    cfg_def.artchan             = cfg.artchan;  % take pairs of electrodes
    cfg_preproc                 = ft_definetrial(cfg_def); 
    %-- remove artefact trial too close from the beginning or end of trial
    cfg_preproc.trl(union(...
        find((cfg_preproc.trl(:,1)-cfg.art_prestim) < 0),...
        find((cfg_preproc.trl(:,2)+cfg.art_poststim) > cfg_preproc.trl(end,2))),:) = [];
    %-- read data using artefact above threshold as trials
    artefact_index{file}        = cfg_preproc.trl;
    cfg_preproc.demean          = 'yes';
    cfg_preproc.baselinewindow  = [cfg.art_prestim -cfg.art_prestim/10]; 
    cfg_preproc.hpfilter        = 'yes';
    cfg_preproc.hpfreq          = .5;
    data_art{file}              = ft_preprocessing(cfg_preproc);
    
    %-- average data
    cfg_artpca                  = [];
    avg_tmp                     = data_art{file};
    avg_tmp.time                = {};
    avg_tmp.trial               = {};
    %-- divide by n trials to reduce memory issue
    for trial = 1:round(length(data_art{file}.trial)/cfg.dividetrial):length(data_art{file}.trial)
        if (trial+round(length(data_art{file}.trial)/cfg.dividetrial)-1) <= length(data_art{file}.trial)
            cfg_artpca.trials   = trial:(trial+round(length(data_art{file}.trial)/cfg.dividetrial)-1);
        else
            cfg_artpca.trials   = trial:length(data_art{file}.trial);
        end
        data_tmp                = ft_timelockanalysis(cfg_artpca,data_art{file});
        avg_tmp.trial{end+1}    = data_tmp.avg;
        avg_tmp.time{end+1}     = data_tmp.time;
    end
    data_art{file}              = ft_timelockanalysis([],avg_tmp);
    
    %-- get continuous signal
    cfg_raw                     = [];
    cfg_raw.continuous          = 'yes';
    cfg_raw.hpfilter            = 'yes';
    cfg_raw.hpfreq              = 1;
    cfg_raw.demean              = 'yes'; % demean on all window
    cfg_raw.baselinewindow      = 'all';
    cfg_raw.dataset             = cfg.dataset{file};
    data_raw{file}              = ft_preprocessing(cfg_raw);
end
%-- append all artefact channels 
if length(cfg.dataset) == 1, 
    data_art                    = data_art{1};
    data_raw                    = data_raw{1};
else
    %-- weighted mean & append
    trial_nb                    = cell2mat(cellfun(@(x) length(x.avg), data_art, 'UniformOutput', false));
    data_art                    = ft_appenddata([], data_art{:});
    data_raw                    = ft_appenddata([], data_raw{:});
    data_art.avg                = [];
    data_raw.concat             = [];
    for file = 1:length(cfg.dataset)
        data_art.avg(:,:,file)  = data_art.trial{file} * (trial_nb(file) / sum(trial_nb));
        data_raw.concat         = cat(2,data_raw.concat,data_raw.trial{file});
    end
    data_art.avg                = squeeze(sum(data_art.avg,3));
    data_art.time               = data_art.time{1};
end

%-- apply pca
data_raw.components             = 0*data_raw.concat;                     % initialize?
data_art.all_comp                    = {};

%-- apply pca independlty for each sensor
for chantype = 1:length(cfg.chantypes)
    disp(['channel type ' num2str(chantype) '/' num2str(length(cfg.chantypes)) ]);
    data_art.all_comp{chantype}      = princomp(data_art.avg(cfg.chantypes{chantype},:)');
    data_raw.components(cfg.chantypes{chantype},:) = data_art.all_comp{chantype}'*data_raw.concat(cfg.chantypes{chantype},:);
    %-- check whether thec components correlate with art channel (to see whether we're capturing some other thing)
    %---- derivative on artchannels
    if cfg.artchan == 1
        deriv_avg               = data_art.avg(cfg.artchan,:)';
        deriv                   = data_raw.concat(cfg.artchan,:)';
    elseif iseven(length(cfg.artchan))
        deriv_avg               = nanmean(...
            data_art.avg(cfg.artchan(1:2:(end-1)),:) - ...
            data_art.avg(cfg.artchan(2:2:end),:))';
        deriv                   = nanmean(...
            data_raw.concat(cfg.artchan(1:2:(end-1)),:) - ...
            data_raw.concat(cfg.artchan(2:2:end),:))';
    else
        error('too many artefact channels or uneven number of channel to apply derivative');
    end
    %---- actual correlation between component space and artchan/derivative of artchan
    disp('computes correlation with art ...');
    data_art.r{chantype}                = abs(corr(data_raw.components(cfg.chantypes{chantype},:)', deriv));
    
    %-- will remove the components that correlates the most with the art
    data_art.rm_components{chantype}    = find(data_art.r{chantype}(1:cfg.max_comp) >= cfg.corr_thresh*nanstd(data_art.r{chantype}));
    %-- and keep the other
    data_art.keep_components{chantype}  = setdiff(1:size(data_art.all_comp{chantype},1),data_art.rm_components{chantype});
    
    data_art.clear_comp{chantype}       = zeros(length(data_art.rm_components{chantype}),size(data_art.all_comp{chantype},1));
    data_art.clear_comp{chantype}       = cat(2,data_art.clear_comp{chantype}',data_art.all_comp{chantype}(:,data_art.keep_components{chantype}))';
    
    
    %-- watch effect on average art
    data_art.avg_comps{chantype}     = ...
        data_art.all_comp{chantype}*...
        data_art.clear_comp{chantype}*...
        data_art.avg(cfg.chantypes{chantype},:);
   
    %-- plotting functions
    if strcmp(cfg.plot,'yes')
        %-- art
        subplot(6,length(cfg.chantypes),length(cfg.chantypes)+chantype);hold on;
        title('art');
        plot(data_art.time,deriv_avg);
        %-- mean artefact
        subplot(3,length(cfg.chantypes),length(cfg.chantypes)+chantype);hold on;
        title(['artifact chans ' num2str(chantype)]);
        plot(data_art.time,data_art.avg(cfg.chantypes{chantype},:),'r'); % average artifact
        plot(data_art.time,data_art.avg_comps{chantype}','g');                % corrected artefact
        legend({'avg artifact', 'corrected artefact'});
        %-- correlation in time between art and principal components
        subplot(3,length(cfg.chantypes) ,2*length(cfg.chantypes)+chantype);hold on;
        title('correlation PC & art')
        scatter(1:length(data_art.r{chantype}), data_art.r{chantype}, '+');
        plot([.5 length(data_art.r{chantype})], repmat(cfg.corr_thresh*nanstd(data_art.r{chantype}),1,2), 'g--');
        plot([cfg.max_comp cfg.max_comp], [0 1], 'g--');
        set(gca,'xscale', 'log');
    end
end

%-- save components
output.chantypes    = cfg.chantypes;
output.all_comp     = data_art.all_comp;
output.clear_comp   = data_art.clear_comp;
output.remove       = data_art.rm_components;
output.trl          = cfg;
output.trial_nb     = trial_nb;
output.artefact_index = artefact_index;

end
