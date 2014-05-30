function [Y cfg] = ft_jr_ecg_ica(cfg,dataset)
% [Y cfg] = ft_jr_ecg_ica(cfg [,dataset])
%--------------------------------------------------------------------------
% Automatic ECG artefact correction based on a combination of ICA and PCA
% % following a series of step:
%    1. Preprocessing
%    2. Remove bad channels
%    3. Normalize channels
%    4. Remove bad trials
%    5. Run fast ica
%    6. Rotate data in C space
%    7. Get nlag derivative
%    8. Identify ECG component
%    9. Construct trials for all chans
%   10. Export to fieldtrip structure
%--------------------------------------------------------------------------
% See below for all default parameters
%--------------------------------------------------------------------------
% Requires following libraries:
%   - Fieldtrip 2011 or more recent
%   - Matlab stats toolbox (functions: trimmean, mad, nanmean, nanstd)
%--------------------------------------------------------------------------
% (c) JeanRémi King 2012, jeanremi.king [at] gmail.com
%--------------------------------------------------------------------------
% 
if nargin == 1 
    if ~isfield(cfg,'dataset'), error('indicate a dataset'); 
    else
        print('0. Read data', 15,55);
        %------ Read continuous data
        if ~isfield(cfg,'read')
            cfg.read.continuous                 = 'yes';
        end
        cfg.read.dataset = cfg.dataset;
        evalc('dataset = ft_preprocessing(cfg.read);');
    end
end
%% default parameters
%------ general
if ~isfield(cfg,'void'),            cfg.void            = true;         end % display feedbacks
global void; void = cfg.void;
if ~isfield(cfg,'plot'),            cfg.plot            = true;         end % display figures
if ~isfield(cfg,'preproc'),         cfg.preproc         = [];           end % preprocessing
if ~isfield(cfg,'ica'),             cfg.ica             = [];           end % independent component analysis
%------ trial definition
if ~isfield(cfg.preproc, 'trialsdef'), 
    cfg.preproc.trialdef.eventtype      = 'STI101';
    cfg.preproc.trialdef.eventvalue     = 1:256;
    cfg.preproc.trialdef.prestim        = 2;
    cfg.preproc.trialdef.poststim       = 2;
end
%------ preprocessing
if ~isfield(cfg.preproc, 'channel'),    cfg.preproc.channel = 'MEG';    end
if ~isfield(cfg.preproc, 'continuous'), cfg.preproc.continuous = 'yes'; end
if ~isfield(cfg.preproc, 'detrend'),    cfg.preproc.detrend = 'yes';    end
if ~isfield(cfg.preproc, 'demean'),     cfg.preproc.detrend = 'demean'; end
%------ resampling
if ~isfield(cfg, 'resample'), 
    cfg.resample.resamplefs             = 150;
    cfg.resample.detrend                = 'no';
end
%------ cleaning
if ~isfield(cfg,'chan_threshold'),      cfg.chan_threshold      = 5;   end
if ~isfield(cfg,'sample_threshold'),    cfg.sample_threshold    = 5;   end

%------ ICA
if ~isfield(cfg.ica,'method'),          cfg.ica.method          = 'fastica';end
if ~isfield(cfg.ica,'numcomponent'),    cfg.ica.numcomponent    = 100;      end
%------ ECG identification
if ~isfield(cfg,'prestim'),             cfg.prestim   = .5;                 end % final trial segmentation
if ~isfield(cfg,'poststim'),            cfg.poststim  = 1;                  end % final segmentation
if ~isfield(cfg,'bsl'),                 cfg.bsl       = [-.500 -.200];      end % baseline correction time
if ~isfield(cfg,'ecg'),
    cfg.ecg                             = [];
    cfg.ecg.sample_threshold            = cfg.sample_threshold;   % artefact threshold, in mad across samples
    cfg.ecg.prestim                     = cfg.prestim;   % artefact trial segmentation
    cfg.ecg.poststim                    = cfg.poststim;    % trial segmentation
    cfg.ecg.bpm                         = 70;   % heart beat per minute
    cfg.ecg.delay                       = 1;    % add one second delay to ensure repeated component
    cfg.ecg.refractory                  = .5;   %prevent too close-in-time components
    cfg.ecg.lag                         = .02;  % n-lag derivative, in ms
    cfg.ecg.timelock_threshold          = 3;    % selected components must be n mad above others components
    cfg.ecg.bsl                         = cfg.bsl; % baseline correction a
end
if ~isfield(cfg,'artprop'),             cfg.artprop = 3/4;                  end % train test proportion of trials
%% Preprocessing: trials, detrend, resample
%-- Define trials based on triggers
print('1. Preprocessing',15,55);
print('1. Define trials',20,50);
cfg.preproc.dataset     = cfg.dataset;
evalc('cfg.preproc      = ft_definetrial(cfg.preproc);');
cfind = @(string,carray) find(cell2mat(cellfun(@(x) ~isempty(strfind(x,string)),carray,'UniformOutput',false)));
% fast manual trigger finding
trigger                 = dataset.trial{1}(cfind(cfg.preproc.trialdef.eventtype,dataset.label),:);
trigger                 = trigger(2:end) - trigger(1:(end-1));
trl                     = find(ismember(trigger,cfg.preproc.trialdef.eventvalue))';
trl(:,1)                = trl(:,1) - (cfg.preproc.trialdef.prestim)*dataset.fsample;
trl(:,2)                = trl(:,1) + (cfg.preproc.trialdef.prestim)*dataset.fsample;
trl(:,3)                = trl(:,2) - trl(:,1);
cfg_tmp                 = [];
cfg_tmp.trl             = trl;
data                    = ft_redefinetrial(cfg_tmp,dataset);
%-- select random trials to get at most 5 minute of data
triggers                = cfg.preproc.trl;
cfg.preproc.trl         = cfg.preproc.trl(randperm(length(cfg.preproc.trl)),:); 
cfg.preproc.trl         = cfg.preproc.trl(1:(5*60/(cfg.preproc.trialdef.poststim+cfg.preproc.trialdef.prestim)),:); 
evalc('data             = ft_preprocessing(cfg.preproc,data);');

%-- resample
print('2. Resample',20,50);
evalc('data             = ft_resampledata(cfg.resample, data);');    % trial based
evalc('dataset          = ft_resampledata(cfg.resample, dataset);'); % continuous
cfg.fsample             = cfg.resample.resamplefs;

%% Remove bad channels
print('2. Remove artefacted channels',15,55);
if ~isfield(cfg,'chantypes'),           cfg.chantypes = 1:length(data.label);          end
chantypes               = cfg.chantypes;
X                       = [data.trial{:}];
badchans                = [];
for c = 1:length(chantypes)
    Xmad                = mad(X(chantypes{c},:),1,2);
    badchans            = cat(2,badchans,chantypes{c}(Xmad>cfg.chan_threshold));
end
cfg.rm_badchan          = badchans;

%% Normalization factor
print('3. Normalize channels',15,55);
n                       = mad(X,1,2);
data.trial              = cellfun(@(x) x./repmat(n,1,size(x,2)), data.trial, 'UniformOutput', false);

%% remove bad trials
print('4. Remove bad trials',15,55);
threshold   = cfg.sample_threshold*mad(reshape([data.trial{:}],1,[]));
values      = cell2mat(cellfun(@(x) std(x(:)), data.trial, 'UniformOutput', false));
data.trial  = data.trial(values<threshold);

%% Run ICA
print('5. Run fast ica',15,55);
chans                   = cellfun(@(x) reshape(x,1,[]), chantypes, 'UniformOutput', false);
chans                   = [chans{:}];
chans(badchans)         = [];
cfg.ica.channel         = setdiff(chans,badchans);
evalc('comp             = ft_componentanalysis(cfg.ica, data);');
% clear data;
% %-- plot components
% cfg_plot = [];
% cfg_plot.component = [1:20];       % specify the component(s) that should be plotted
% cfg_plot.layout = 'neuromag306all.lay';
% cfg_plot.comment   = 'no';
% ft_databrowser(cfg_plot,comp);

%% Rotation
print('6. Rotate data in C space',15,55);
%-- filter
cfg_tmp                 = [];
cfg_tmp.detrend         = 'yes';
cfg_tmp.demean          = 'yes';
cfg_tmp.hpfilter        = 'yes';
cfg_tmp.hpfreq          = .8;
cfg_tmp.channel         = cfg.ica.channel;
evalc('X                = ft_preprocessing(cfg_tmp, dataset);');
%-- select channels
X                       = X.trial{1};
%-- unnormalize components
Cweights                = comp.topo .* repmat(n,1,size(comp.topo,2));
%-- rotate in component space
Xc                      = norm(Cweights)'*(norm(X));

%% Component transform
%-- transform into derivative
print('7. Get nlag derivative',15,55);
lag                     = ceil(cfg.ecg.lag*cfg.resample.resamplefs);
Xc                      = abs(Xc(:,lag:end)-Xc(:,1:(end-lag+1)));

%% ECG identification
print('8. Identify ECG component',15,55);
fsample                 = cfg.resample.resamplefs;
refractory              = cfg.ecg.refractory;                   % minimum time separating two tirals
bpm                     = cfg.ecg.bpm;
n_comps                 = size(Xc,1);
delay                   = cfg.ecg.delay;
pre                     = delay-cfg.ecg.prestim;        % time taken before artefact threshold:
post                    = delay+cfg.ecg.poststim;       % time taken after artefact threshold
t0                      = floor(pre*fsample);
n_samples               = ceil(1+(post-pre)*fsample);
bsl                     = delay+cfg.ecg.bsl;
bsl                     = 1-t0+((floor(bsl(1)*fsample)):((ceil(bsl(2)*fsample))));% time for baseline correction
expected_freq           = 1-2*bpm/60/fsample; % frequency of occurence across samples, time 2 samples above threshold
[Cmed Cmad]             = deal(NaN(n_comps,n_samples));%initialize
[all_thresholds all_trl]= deal(cell(n_comps,1));
print('1. Components trials alignment',20,50);
for c = 1:n_comps
    bar(c,n_comps);
    %----- clean components
    threshold           = nanmedian(abs(Xc(c,:)))+cfg.ecg.sample_threshold*mad(abs(Xc(c,:)));
    good_samples        = abs(Xc(c,:))<threshold;
    
    %------ find threshold to get between ~1 threshold per second
    [h bins]            = hist(abs(Xc(c,good_samples)),1000);
    
    %------ cumulative probability of being higher than a given threshold;
    p                   = cumsum(h,2)./repmat(sum(h,2),1,1000);
    
    %------ get corresponding threshold for each component
    threshold           = bins(find(p>=expected_freq,1));
    all_thresholds{c}   = threshold;
    
    %------ align trials
    trl                 = find(abs(Xc(c,:))>=threshold);
    
    %------ save all trl for later use
    all_trl{c}          = trl;
    %-- remove trials Yside recording
    bad_trials          = [find(floor(trl+pre*fsample)<1) find(ceil(trl+post*fsample)>size(Xc,2))];
    trl(bad_trials)     = [];
    %-- remove consecutive trials
    bad_trials          = (trl(2:end)-trl(1:(end-1)))<(refractory*fsample);
    trl(bad_trials)     = [];
    
    %-- if found at least 25% of the component that should have occured
    if length(trl)>(.25*size(Xc,2)/(bpm/60*fsample))
        trials= NaN(length(trl),n_samples);
        for t = 1:length(trl)
            %-- select trial
            t_sel       = floor(trl(t)+pre*fsample):ceil(trl(t)+post*fsample);
            trials(t,:) = Xc(c,t_sel)- nanmean(Xc(c,bsl));
        end
        %---- save trials
        Cmed(c,:)       = norm(trimmean(trials,90,'round',1));
        Cmad(c,:)       = norm(mad(trials,90)); % mad indicates subsequent
        % events (positive mad) or highly synchronized events (negative mad).
        % Abs norm mad hence identify both phenomenon.
    else
        Cmed(c,:)       = 0;
        Cmad(c,:)       = 0;
    end
end
print('2. Gest best component',20,50);
%-- get best component
Cmedmad                 = Cmad.*Cmed;
[scores order]          = sort(max(Cmedmad,[],2),'descend');
bad_c                   = unique([find(isnan(scores)) find(isinf(scores))]);
order                   = order([setdiff(1:length(order), bad_c) bad_c]);
scores                  = scores([setdiff(1:length(order), bad_c) bad_c]);
if scores(1)<(cfg.ecg.timelock_threshold*nanmedian(scores))
    print('Warning: did not converge towards a single component!',25, 45);
end
%---- keep selected component
artchan                 = Xc(order(1),:);
clear Xc
trl                     = all_trl{order(1)};

%% Get trial data for all channels
print('9. Construct trials for all chans',15,55);
%-- parameters
pre                     = -cfg.prestim;    % time taken before artefact threshold:
post                    = cfg.poststim;   % time taken after artefact threshold
t0                      = floor(pre*fsample);
time                    = (pre*fsample):(post*fsample);
n_samples               = 1+ceil((post-pre)*fsample);
bsl                     = cfg.bsl;
bsl                     = 1-t0+((floor(bsl(1)*fsample)):((ceil(bsl(2)*fsample))));% time for baseline correction

%-- remove trials Yside recording
bad_trials              = [find(floor((trl+pre*fsample))<1) find(ceil(trl+post*fsample)>size(X,2))];
trl(bad_trials)         = [];

%-- remove consecutive trials
bad_trials              = (trl(2:end)-trl(1:(end-1)))<((post-pre)*fsample/3);
trl(bad_trials)         = [];

%-- align trials
trials                  = NaN(size(X,1),n_samples,length(trl));
artchan_trl             = NaN(1,n_samples,length(trl));
for t = 1:length(trl)
    bar(t,length(trl));
    t_sel               = floor(trl(t)+pre*fsample):ceil(trl(t)+post*fsample);
    %-- select trial
    trial               = X(:,t_sel);
    %-- apply baseline
    trials(:,:,t)       = trial - repmat(nanmean(trial(:,bsl),2),1,length(t_sel));
    %-- 
    artchan_trl(:,:,t)  = artchan(:,t_sel);
end


%%  Output
%-- transform trials data into fieldtriplike structure:
print('10. Export ',15,55);
Y.Cweights              = Cweights;
Y.Corder                = order;
Y.Cscores               = scores;
Y.Cmedmad               = Cmedmad;
Y.Cthresholds           = all_thresholds;
Y.Ctrl                  = all_trl{c};
Y.data_artchan          = artchan;
Y.artchan_trl           = artchan_trl;
Y.time                  = time;
Y.trial_info            = [trl-pre trl+post (post-pre)]*fsample;
for n = 1:size(trials,1)
    bar(t,size(trials,1));
    Y.trial{n}          = squeeze(trials(:,:,n));
end
rnd_trl                  = randperm(n);
Y.trl_sel.train          = rnd_trl(1:round(cfg.artprop*n));
Y.trl_sel.test           = rnd_trl((n-round((1-cfg.artprop)*n)):n);
Y.all_trl               = trl;
Y.cfg                   = cfg;

%-- randomize trials and create test and train sets: note that these won't
%be fully independent as the ICA was applied on the entire dataset.
end

function print(txt,left,right)
%-- automatic printing
if nargin == 1,
    % center
    left = 72/2 - floor(length(txt)/2) -1;
    right = 72/2 + ceil(length(txt)/2) -1;
end
global void
if void
    disp([repmat('-',1,left) ' ' txt ' ' repmat('-',1,right-length(txt))]);
end
end

function bar(x,n)
%progressing bar
global void
if void
    bar_ = {'','*', '.\n'};
    fprintf(bar_{1+double((ceil(72*(x+1)/n)-ceil(72*x/n))==1)+double(x==n)+(n<72)*(x==n)});
end
end

function y = norm(x)
y= (x-trimmean(x(:),90))./mad(x(:));
end
