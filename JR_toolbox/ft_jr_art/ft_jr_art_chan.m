function [data_out cfg]= ft_jr_art_chan(cfg)

%% Default parameters
if ~isfield(cfg,'artchan'),     error('cfg.artchan is needed!');            end
if ~isfield(cfg,'void'),        cfg.void        = true;                     end % void
if ~isfield(cfg,'threshold'),   cfg.threshold   = 3.5;                      end % threshold to define trial based on art channel
if ~isfield(cfg,'dataset'),     error('needs cfg.dataset');                 end
if ~isfield(cfg,'artchan'),     error('cfg.artchan is needed!');            end
if ~isfield(cfg,'chantypes'),   cfg.chantypes   = {1:hdr.nChans};           end % cells dividing types of sensors (gradiometers, etc)
if ~isfield(cfg,'artprop'),     cfg.artprop     = 3/4;                      end % max number of artefacts to keep
if ~isfield(cfg,'memory'),      cfg.memory      = 'low';                    end % memory needs
if ~isfield(cfg,'badchan_threshold'),cfg.badchan_threshold = 4;             end % automatic bad channels removal in n * MAD
if ~isfield(cfg,'preproc'),     cfg.preproc = [];                           end % default preproc, see below
cfg.trialfun                    = 'ft_jr_art_trialfun';                         % get heart beat as triggers

%% Define generic functions
%------ compatibility
% if ~exist('nanmean','file'),    nanmean= @(x,dim) mean(x,dim); warning('Could not find nanmean'); end
% if ~exist('nanstd','file'),     nanstd = @(x,flag,dim) std(x,flag,dim); warning('Could not find nanstd'); end
% if ~exist('trimmean','file'),   trimmean = @(x,percent,flag,dim) nanmean(x,dim); end
% if ~exist('mad','file'),        mad = @(x,flag,dim) nanstd(x,flag,dim); end
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,15) ' ' x ' ' repmat('-',1,55-length(x))]);
else
    print = @(x) false;
end
%------ progressing bar
if cfg.void
    bar_ = {'','*', '.\n'};
    bar = @(x,n) fprintf(bar_{1+mod(ceil(70*(x+1)/n)-ceil(70*x/n),2)+double(x==n)});
else
    bar = @(x,n) false;
end
%% Main

global threshold data_artchan% get threshold back from ft_jr_trialfun
threshold = cfg.threshold;


%% 1. find artefacted moments
print('1. Find artifact timing');
evalc('cfg_art              = ft_definetrial(cfg);') % silence output
if isnan(threshold), data_art = []; warning('could not find artefact'); return; end
cfg.threshold               = threshold; % update threshold 
%-- randomize trials
cfg_preproc                 = cfg.preproc;
cfg_preproc.dataset         = cfg.dataset;
cfg_preproc.trl             = cfg_art.trl(randperm(size(cfg_art.trl,1)),:);


%% 2. Get artefacted ERP
print('2. Preprocess data ');
%-- preproc for all channels
cfg_preproc.continuous = 'no';
if ~isfield(cfg_preproc,'bpfilter'),        cfg_preproc.ppfilter = 'yes';   end
if ~isfield(cfg_preproc,'bpfreq'),          cfg_preproc.ppfreq = [1 20];    end
if ~isfield(cfg_preproc,'demean'),          cfg_preproc.demean = 'yes';     end
if ~isfield(cfg_preproc,'baselinewindow'),  cfg_preproc.baselinewindow= [-.200 0]; end
evalc('data_art             = ft_preprocessing(cfg_preproc);'); % silence output
cfg.fsample                 = data_art.fsample;
% select independent trials
n                           = size(cfg_art.trl,1);
data_art.trl_sel.train      = 1:round(cfg.artprop*n);
data_art.trl_sel.test       = (n-round((1-cfg.artprop)*n)):n;

%-- sample for component
print('3. Preprocess component ');
cfg_tmp             = [];
cfg_tmp.time        = (1:length(data_artchan))/data_art.fsample;
cfg_tmp.label       = {'artchan'};
cfg_tmp.fsample     = data_art.fsample;
data_artchan_trl    = data2ft(cfg_tmp,reshape(data_artchan,1,[]));
evalc('data_artchan_trl     = ft_preprocessing(rmfield(cfg_preproc,''dataset''),data_artchan_trl);'); 
evalc('data_artchan_trl     = ft_redefinetrial(rmfield(cfg_preproc,''dataset''),data_artchan_trl);'); 
data_art.time               = data_art.time{1};



%% 3. Automatic bad channel detection
print('4. Automatic bad channel detection');
if ~isfield(cfg,'rm_badchan')
    switch cfg.memory
        case 'low'
            %-- get median absolute deviation for each channel in order to
            %prevent crap channel from exploding the variance
            chan_mad = NaN(length(data_art.label),1);
            n_chans = length(data_art.label);
            for channel = 1:n_chans
                bar(channel,n_chans);
                chan_mad(channel)    = mad(cell2mat(cellfun(@(x) x(channel,:), data_art.trial(data_art.trl_sel.train), 'UniformOutput', false)),1,2);
            end
            
            %-- find outlier channel
            n_chantypes = length(cfg.chantypes);
            rm_badchan = cell(n_chantypes);
            for chantype = 1:n_chantypes
                rm_badchan{chantype} = intersect(...
                    cfg.chantypes{chantype},...
                    find(chan_mad>(cfg.badchan_threshold*median(chan_mad))));
            end
            rm_badchan = [rm_badchan{:}];
            cfg.rm_badchan  = rm_badchan;
        otherwise, error('unknown cfg.memory option');
    end
end


%% Output
data_out            = [];
data_out.time       = data_art.time;
data_out.trial      = data_art.trial;
data_out.data_artchan=data_artchan;
data_out.artchan_trl= reshape(cell2mat(data_artchan_trl.trial),[1 length(data_artchan_trl.trial{1}) length(data_artchan_trl.trial)]);
data_out.all_trl    = cfg_preproc.trl;
data_out.trl_sel    = data_art.trl_sel;
cfg.preproc         = cfg_preproc;
cfg.sample_threshold= threshold;

clear data_artchan threshold;