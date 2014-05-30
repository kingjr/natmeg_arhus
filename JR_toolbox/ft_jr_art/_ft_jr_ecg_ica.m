function [Y cfg] = ft_jr_ecg_ica(cfg,dataset)
% [Y cfg] = ft_jr_ecg_ica(cfg [,dataset])
%--------------------------------------------------------------------------
% Automatic ECG artefact correction based on a combination of ICA and PCA
%--------------------------------------------------------------------------
%   Input parameter           default value
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       cfg.void                            Display feedbacks
%       cfg.plot                            Display figures
%       cfg.bpfilter           'yes'
%       cfg.bpfreq              [1 20]
%
%       cfg.chan_trigger:       'STI101'
%       cfg.chan_meg:           'MEG'
%       cfg.chan_ecg:           'ECG'
%       cfg.rm_art_chan_th      10          Threshold, in mad, above
%                                           which channel will be excluded
%                                           from ICA
%       cfg.rm_art_sample_th    4           Idem as above for samples
%       cfg.rm_art_pad          20          Removed n neighboring samples
%       cfg.maxsample       	10000       Number of samples used in ICA
%       cfg.chantypes           {grad mag}  Specify distinct types of
%                                           chans
%       cfg.norm_chan           10.^[11 12] Avoid numerical errors in ICA
%       cfg.sample_threshold    10           threshold, in mad, above wich
%                                           abs(data_comp) will be set down
%                                           to 0;
%       cfg.bpm                 70          Ecg beats per minute
%       cfg.refractory          .5          ecg refractory period
%       cfg.poststim            1.000       Timing around artefact findings
%                                           in s (takes 1 second after!)
%       cfg.prestim             .500        Time before threshold
%       cfg.bsl                 -[.500 .700]Baseline time for bsl
%                                           correction on redefinition of
%                                           trials, in s
%       cfg.timelock_threshold  3           Threshold to consider component
%                                           likely reflect single ECG
%                                           components, in mad.
%       cfg.realign_threshold   .010        max time difference between
%                                           identical events founds across
%                                           ICAs, in s
%       cfg.artprop             3/4         train test trials proportion
%       [unused: cfg.n_trials    200        Maximum number of trials to be
%                                           considered]
%
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       [dataset]:              dataset as a fieldtripstructure. If not
%                               provided, cfg.dataset must be given.
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Output:
%       Y
%       cfg
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%--------------------------------------------------------------------------
% Requires following libraries:
%   - Fieldtrip 2011 or more recent
%   - Matlab stats toolbox (functions: trimmean, mad, nanmean, nanstd)
%--------------------------------------------------------------------------
% (c) JeanRémi King 2012, jeanremi.king [at] gmail.com
%--------------------------------------------------------------------------
% 
if nargin == 1 && ~isfield(cfg,'dataset'), error('indicate a dataset'); end
%% default parameters
if ~isfield(cfg,'void'),            cfg.void            = true;         end % display feedbacks
if ~isfield(cfg,'plot'),            cfg.plot            = true;         end % display figures

if ~isfield(cfg,'chan_trigger'),    cfg.chan_trigger    = 'STI101';     end % trigger channel
if ~isfield(cfg,'chan_meg'),        cfg.chan_meg        = 'MEG';        end % MEG channels
if ~isfield(cfg,'chan_ecg'),        cfg.chan_ecg        = 'ECG';        end % ECG channel
if ~isfield(cfg,'bpfilter'),        cfg.bpfilter        = 'yes';        end % ECG channel
if ~isfield(cfg,'bpfreq'),          cfg.bpfreq          = [1 20];       end % ECG channel

if ~isfield(cfg,'rm_art_chan_th'),  cfg.rm_art_chan_th  = 15;           end % threshold, in mad, above which channel will be excluded from ICA
if ~isfield(cfg,'rm_art_sample_th'),cfg.rm_art_sample_th= 4;            end % threshold, in mad, above which samples will be excluded from ICA
if ~isfield(cfg,'rm_art_pad'),      cfg.rm_art_pad      = 20;           end % padding, in sample, to remove time to be excluded from ICA based on threshold

if ~isfield(cfg,'maxsample'),       cfg.maxsample       = 100000;        end % max number of samples used in ICA;

if ~isfield(cfg,'chantypes'),       cfg.chantypes       = {[1:3:306 2:3:306]; 3:3:306}; end % gradiometers and magnetometers
if ~isfield(cfg,'norm_chan'),       cfg.norm_chan       = 10.^[11 12];  end % normalization factor for each chantype

if ~isfield(cfg,'sample_threshold'),cfg.sample_threshold= 10;            end % threshold, in mad, above wich abs(data_comp) will be set down to 0;
if ~isfield(cfg,'bpm'),             cfg.bpm             = 70;           end % ecg beats per minute
if ~isfield(cfg,'refractory'),      cfg.refractory      = .5;           end % ecg refractory period

if ~isfield(cfg,'poststim'),        cfg.poststim        = 1.500;        end % timing around artefact findings, in s
if ~isfield(cfg,'prestim'),         cfg.prestim         = .700;         end %
if ~isfield(cfg,'bsl'),             cfg.bsl             = -[.500 .200]; end % baseline time for bsl correction on redefinition of trials, in s
if ~isfield(cfg,'n_trials'),        cfg.n_trials        = 200;          end % maximum number of trials to be considered
if ~isfield(cfg,'timelock_threshold'),cfg.timelock_threshold = 3;       end % threshold to consider component likely reflect single ECG components as compared to other components, in mad.

if ~isfield(cfg,'realign_threshold'),cfg.realign_threshold = .010;      end % max time difference between identical events founds across ICAs, in s
if ~isfield(cfg,'artprop'),         cfg.artprop     = 3/4;              end % max number of artefacts to keep


%% Define generic functions
global void;void=cfg.void;%for bar and print
if ~exist('fastica.m'),
    pwd = which('fieldtripdefs.m');
    addpath([pwd(1:(length(pwd)-length('compat/fieldtripdefs.m'))) '/external/fastica/']);
end
% select trials with specific sensors from fieldtrip structure:
sel_sensors = @(data,chan) reshape(cell2mat(cellfun(@(x) x(chan,:), data, 'UniformOutput',false)),[1 length(data{1}) length(data)]);

%% read data
print('1. Read data',15,55);
if nargin == 1
    cfg_dataset         = [];
    cfg_dataset.dataset = cfg.dataset;
    evalc('dataset = ft_preprocessing(cfg_dataset);');
end
cfg.fsample = dataset.fsample;

%% select time between first and last trigger
print('2. Define Chan Of Interest & Time o.I.',15,55);
chan_meg    = cfind(cfg.chan_meg,       dataset.label);
chan_sti101 = cfind(cfg.chan_trigger,   dataset.label);
data_sti101 = sel_sensors(dataset.trial, chan_sti101);
start       = find(round(data_sti101)>0,1);
stop        = find(round(data_sti101)>0,1,'last');
time        = start:min(stop-start,round((stop-start)/(length(chan_meg)*100))):stop;

%% MEG filtering
print('3. Filter MEG',15,55);
cfg_bp         = [];
cfg_bp.bpfilter= cfg.bpfilter;
cfg_bp.bpfreq  = cfg.bpfreq;
evalc('dataset = ft_preprocessing(cfg_bp,dataset);');


% select time and channels
for chantype = 1:length(cfg.chantypes)
    print(['[Chantype ' num2str(chantype) ' ...'],15,55);
    chans               = cfg.chantypes{chantype};
    
    print('4. Remove artefacted chans & time',15,55);
    Xc                  = dataset.trial{1}(chans,time);
    [Xc rm_badchan]    = select_clean_data(Xc,...     % continuous data
        cfg.norm_chan(chantype),...                 % normalization factor to avoid numerical errors
        cfg.rm_art_chan_th,...                      % channel threshold removal in mad
        cfg.rm_art_sample_th,...                    % sample threshold removal in mad
        cfg.rm_art_pad);
    Y.rm_badchan{chantype}      = chans(rm_badchan);
    
    %% compute ICA
    print('5. Compute fast ICA',15,55);
    cfg.samples                 = randperm(min(size(Xc,2),cfg.maxsample));
    [unused unused weights]     = fastica(Xc(:,cfg.samples),'maxNumIterations', 500,'verbose','on');
    unused = [];clear unused;
    
    %-- resample
    cfg_resample            = [];
    cfg_resample.resamplefs = 150;
    cfg_resample.detrend    = 'no';
    data_resample            = ft_resampledata(cfg_resample, dataset);
    %-- filter out low frequency
    cfg_preproc = [];
    cfg_preproc.continuous = 'yes';
    cfg_preproc.hpfilter = 'yes';
    cfg_preproc.hpfreq = .5;
    data_preproc = ft_preprocessing(cfg_preproc,data_resample);
    %-- remove bad channels
    X           = data_preproc.trial{1}(chan_meg,:)';
    S           = sum(abs(X));
    threshold   = trimmean(S,90)+cfg.rm_art_chan_th*mad(S);
    sel_chans   = find(S<=threshold);
    rm_badchans = setdiff(1:length(S),sel_chans);
    
    %-- normalize data
    for c = 1:length(chan_meg)
        X(c,:) = norm(X(c,:));
    end
    %------ identify bad samples
    S           = sum(abs(X'));
    threshold   = trimmean(S,90)+cfg.rm_art_sample_th*mad(S);
    %------ identify bad channels
    sel_samples   = find(S<=threshold);
    rm_samples = setdiff(1:length(S),sel_chans);
    %------ output in fieldtrip
    data_preproc.time{1}    = data_preproc.time{1}(sel_samples);
    data_preproc.label      = data_preproc.label(chan_meg(sel_chans));
    data_preproc.trial{1}   = data_preproc.trial{1}(chan_meg(sel_chans),sel_samples);
    cfg_ica                 = [];
    cfg_ica.method          = 'runica';
    comp                    = ft_componentanalysis(cfg_ica, dataset);
    
    %% transform data into component space
    print('6. Rotate data to components space',15,55);
    data_comp                   = weights*dataset.trial{1}(chans,:);
    R = corr(dataset.trial{1}(cfind('ECG',dataset.label),:)',data_comp');
    scatter(1:length(R),R);
    
    %% identify ECG component based on dynamics of averaged signals
    print('7. Identify ECG component',15,55);
    [Cmedmad trl order scores]  = mean_ComponentSpace(data_comp, cfg);
    
    %% Output
    Y.data_artchan{chantype}    = data_comp(order(1),:);
    Y.Cweights{chantype}        = weights(order,:);
    Y.Cmedmad{chantype}         = Cmedmad(order(1),:);
    Y.Cscores{chantype}         = scores;
    Y.Ctrl{chantype}            = trl{order(1)};
    
    print(['... // chantype ' num2str(chantype) ' ]'],15,55);
end


% identify common events across each chantype's ICA
print('8. Redefine trials with components combination',15,55);
%---- concatenate components by multiplying their abs normalized value
data_artchan                    = reshape(cell2mat(...
    cellfun(@(x) norm(x), Y.data_artchan,'UniformOutput',false)),...
    [],length(Y.data_artchan))';
data_artchan                    = abs(prod(data_artchan,1));
Y.data_artchan                  = data_artchan;
%----- redefine trials
threshold                       = nanmedian(data_artchan)+cfg.sample_threshold*mad(data_artchan);
trl                             = find(data_artchan>=threshold);

%----- remove trials Yside recording
bad_trials                      = [...
    find(floor(trl-cfg.prestim*cfg.fsample)<1) ...
    find(ceil(trl+cfg.poststim*cfg.fsample)>size(Y.data_artchan,2))];
trl(bad_trials)                 = [];
%----- remove consecutive trials
bad_trials                      = (trl(2:end)-trl(1:(end-1)))<(cfg.refractory*cfg.fsample);
trl(bad_trials)                 = [];
% clf;plot(data_artchan);hold on;scatter(trl,threshold*ones(length(trl),1),'r')
artchan_trl                     = mean_channelSpace(Y.data_artchan,trl,cfg);
Y.artchan_trl                   = reshape(cell2mat(artchan_trl.trial),[1,length(artchan_trl.trial{1}) length(artchan_trl.trial)]);


%% Get trial data for all channels
print('9. Construct trials for all chans',15,55);
data_trl                        = mean_channelSpace(squeeze(dataset.trial{1}),trl,cfg);

%%  Output
cfg.rm_badchan                  = [Y.rm_badchan{:}];
Y.all_trl                       = data_trl.trial_info;
Y.trial                         = data_trl.trial;

Y.time                          = data_trl.time{1};
%-- randomize trials and create test and train sets: note that these won't
%be fully independent as the ICA was applied on the entire dataset.
n                               = length(data_trl.trial);
rnd_trl                         = randperm(n);
Y.trl_sel.train                 = rnd_trl(1:round(cfg.artprop*n));
Y.trl_sel.test                  = rnd_trl((n-round((1-cfg.artprop)*n)):n);

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

function index = cfind(string,carray)
% index = cfind(string,carray)
% finds strings in a cell array
% (c) JeanRemi King, 2012.
% jeanremi.king [at] gmail [dot] com
index = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,string)),carray,'UniformOutput',false)));
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

function [Xc rm_badchans] = select_clean_data(X,normalization_factor,threshold_chans,threshold_samples,pad)


%% normalization of grad and mag
print('1. Normalize data',20,50)
X           = X.*normalization_factor;

%% set bad channels to 0 for the following ICA
print('2. Identify bad channels',20,50)
%------ identify threshold
S           = sum(abs(X)');
threshold   = trimmean(S,90)+threshold_chans*mad(S);
%------ identify bad channels
sel_chans   = find(S<=threshold);
rm_badchans = setdiff(1:length(S),sel_chans);
if ~isempty(rm_badchans)
    print(['/!\ ' num2str(length(rm_badchans)) ' channel(S) removed from ICA!'],25,45);
    pause(3);
end

% clf;hold on;plot(S);plot([1 length(S)], repmat(threshold,1,2),'r');

%% remove bad times from ICA
print('3. Identify bad timings',20,50)
%------ identify threshold
S           = sum(X(sel_chans,:).^2);
threshold   = trimmean(S,90)+threshold_samples*mad(S);
% clf;hold on;plot(S);plot([1 length(S)], repmat(threshold,1,2),'r');
%------ remove bad timings
Xc          = X;
Xc(setdiff(1:size(Xc,1),sel_chans),:) = 0;
art_time    = find(S>threshold);
bad_time    = unique(reshape(repmat(art_time,pad,1) + repmat(round(linspace(-pad/2,pad/2,pad))',1,length(art_time)),1,[]));
bad_time([find(bad_time<1) find(bad_time>size(Xc,2))]) = [];
Xc(:, bad_time) = [];

%% demean
print('4. Demean',20,50)
Xc              = Xc - repmat(trimmean(Xc,90,'round',2),1,size(Xc,2));

end

function [C all_trl order scores] = mean_ComponentSpace(data_comp, cfg)

%-- parameters
delay           = 1;                    % add one second delay to ensure repeated component
n_comps         = size(data_comp,1);
fsample         = cfg.fsample;
pre             = delay-cfg.prestim;    % time taken before artefact threshold:
post            = delay+cfg.poststim;   % time taken after artefact threshold
t0              = floor(pre*fsample);
n_samples       = ceil(1+(post-pre)*fsample);
bsl             = delay+cfg.bsl;
bsl             = 1-t0+((floor(bsl(1)*fsample)):((ceil(bsl(2)*fsample))));% time for baseline correction
% n_trials        = cfg.n_trials;         % maximum number of trials
bpm             = cfg.bpm;              % ecg beat per min

%----- clean components
print('1. Clean components',20,50);
good_samples = cell(n_comps,1);
for c = 1:n_comps
    bar(c,n_comps);
    threshold = nanmedian(abs(data_comp(c,:)))+cfg.sample_threshold*mad(abs(data_comp(c,:)));
    good_samples{c} = abs(data_comp(c,:))<threshold;
    %     subplot(ceil(sqrt(n_comps)),ceil(sqrt(n_comps)),c);
    %     plot(abs(data_comp(c,:)));hold on;
    %     plot(xlim,[threshold threshold],'r');
    %     axis off
end


%------ find threshold to get between ~1 threshold per second
print(['2. Find thresholds for bpm=' num2str(bpm)],20,50);
[h bins]        = deal(NaN(n_comps,1000));
for c = 1:n_comps
    bar(c,n_comps);
    [h(c,:) bins(c,:)] = hist(abs(data_comp(c,good_samples{c})),1000);
end
expected_freq   = 1-10*bpm/60/fsample; % frequency of occurence across samples, time 10 samples above threshold

%------ cumulative probability of being higher than a given threshold;
p               = cumsum(h,2)./repmat(sum(h,2),1,1000);

%------ get corresponding threshold for each component
thresholds = NaN(n_comps,1);
for c = 1:n_comps
    thresholds(c) = bins(c,find(p(c,:)>=expected_freq,1));
%     clf;plot(abs(data_comp(c,:)));hold on;plot(xlim,[thresholds(c) thresholds(c)],'r');
end

%% define trials
print('3. Align beats to trial',20,50);
[Cmed Cmad]     = deal(NaN(n_comps,n_samples));%initialize
all_trl         = cell(n_comps,1);
for c = 1:n_comps
    bar(c,n_comps);
    trl         = find(abs(data_comp(c,:))>=thresholds(c));
    all_trl{c}  = trl;
    %-- remove trials Yside recording
    bad_trials  = [find(floor(trl+pre*fsample)<1) find(ceil(trl+post*fsample)>size(data_comp,2))];
    trl(bad_trials) = [];
    %-- remove consecutive trials
    bad_trials  = (trl(2:end)-trl(1:(end-1)))<(cfg.refractory*fsample);
    trl(bad_trials) = [];
    %     %-- randomize trials
    trl         = trl(randperm(length(trl)));
    %-- align trials
    %trl         = trl(1:min(n_trials,length(trl)));

    %-- if found at least 25% of the component that should have occured
    if length(trl)>(.25*size(data_comp,2)/(bpm/60*fsample))
        trials = NaN(length(trl),n_samples);
        for t = 1:length(trl)
            %-- select trial
            t_sel           = floor(trl(t)+pre*fsample):ceil(trl(t)+post*fsample);
            trials(t,:) = data_comp(c,t_sel);
            %-- apply baseline
            trials(t,:) = trials(t,:) - nanmean(trials(t,bsl));
        end

        Cmed(c,:) = norm(median(abs(trials)));
        Cmad(c,:) = abs(norm(mad(trials,90))); % mad indicates subsequent
        % events (positive mad) or highly synchronized events (negative mad).
        % Abs norm mad hence identify both phenomenon.
        %         imagesc(trials);pause;
    else
        Cmed(c,:) = 0;
        Cmad(c,:) = 0;
    end
end
% imagesc(Cmed);
% output trl in ft format

%-- get trials in channel space
print('4. Identify most likely time-locked component',20,50);
C=Cmad.*Cmed;
[scores order] = sort(max(C,[],2),'descend');

bad_c   = unique([find(isnan(scores)) find(isinf(scores))]);
order   = order([setdiff(1:length(order), bad_c) bad_c]);
scores  = scores([setdiff(1:length(order), bad_c) bad_c]);
% subplot(221);imagesc(Cmed);
% subplot(222);imagesc(Cmad);
% subplot(223);imagesc(Cmed.*Cmad);
% subplot(224);scatter(order,scores,[],[1 0 0;zeros(length(order)-1,3)], 'filled');

if scores(1)<(cfg.timelock_threshold*nanmedian(scores))
    print('Warning: did not converge towards a single component!',25, 45);
    pause(1);
end

end

function ft = mean_channelSpace(data,trl,cfg)

%-- parameters
fsample         = cfg.fsample;
pre             = -cfg.prestim;    % time taken before artefact threshold:
post            = cfg.poststim;   % time taken after artefact threshold
t0              = floor(pre*fsample);
n_samples       = 1+ceil((post-pre)*fsample);
bsl             = cfg.bsl;
bsl             = 1-t0+((floor(bsl(1)*fsample)):((ceil(bsl(2)*fsample))));% time for baseline correction
% n_trials        = cfg.n_trials;         % maximum number of trials

%-- remove trials Yside recording
print('1. Build trials',20,50);
bad_trials          = [find(floor((trl-pre*fsample))<1) find(ceil(trl+post*fsample)>size(data,2))];
trl(bad_trials)     = [];
%-- remove consecutive trials
bad_trials          = (trl(2:end)-trl(1:(end-1)))<((post-pre)*fsample)/3;
trl(bad_trials)     = [];
% %-- randomize trials
% trl                 = trl(randperm(length(trl)));
% %-- align trials
% trl                 = trl(1:min(n_trials,length(trl)));
trials              = NaN(length(trl),size(data,1),n_samples);
for t = 1:length(trl)
    bar(t,length(trl));
    t_sel           = floor(trl(t)+pre*fsample):ceil(trl(t)+post*fsample);
    %-- select trial
    trial           = data(:,t_sel);
    %-- apply baseline
    trials(t,:,:)   = trial - repmat(nanmean(trial(:,bsl),2),1,length(t_sel));
end

time = (pre*cfg.fsample):(post*cfg.fsample); % time

%-- transform trials data into fieldtriplike structure:
print('2. Export to fieldtrip structure',20,50);
for t = 1:size(trials,1)
    bar(t,size(trials,1));
    ft.trial{t}             = squeeze(trials(t,:,:));
    ft.time{t}              = time;
end
ft.trial_info = [trl-pre*fsample trl+post*fsample (post-pre)*fsample];
end
