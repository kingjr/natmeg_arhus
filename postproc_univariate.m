function [data_erf data_p] = postproc_univariate(data,y,cfg)
if ~exist('cfg','var'), cfg = []; end
if ~isfield(cfg, 'test'),   cfg.test = 'ranksum'; end
if ~isfield(cfg, 'out'),   cfg.out = 'fieldtrip'; end

% generic function
get_trials = @(x) reshape(cell2mat(x.trial),[size(x.trial{1}) length(x.trial)]); % fast way of getting ft_trials
selchan = @(c,s) find(cell2mat(cellfun(@(x) ~isempty(strfind(x,s)),c,'uniformoutput', false))==1);


time    = data.time{1};
dat     = permute(get_trials(data),[3 1 2]);

% channel selection
%if ~isfield(cfg, 'channel'),    cfg.channel = selchan(data.label,'MEG'); end
if ~isfield(cfg, 'channel'),    cfg.channel = 1:length(data.label); end

% time of interest
%if ~isfield(cfg, 'toi'),    cfg.toi = find(time>-.100 & time < .600); end
if ~isfield(cfg, 'toi'),    cfg.toi = 1:length(data.time{1}); end


% face 1 versus face 2
switch cfg.test
    case 'ranksum'
        % mann u whitney test, which is the equivalence of the independent
        % sample t-test for non parametric statistics
        p = ranksum_fast(dat(y==1,cfg.channel,cfg.toi),dat(y==2,cfg.channel,cfg.toi));
    case 'ttest2'
        [~, p] = ttest2(dat(y==1,cfg.channelchannel,cfg.toi),dat(y==2,cfg.channel,cfg.toi));
end

% compute mean
erf = squeeze(...
    trimmean(dat(y==1,cfg.channel,cfg.toi),90,'round',1) - ...
    trimmean(dat(y==2,cfg.channel,cfg.toi),90,'round',1));

switch cfg.out
    case 'fieldtrip'
        % transform into fieldtrip strcture
        % p value
        data_p = data;
        data_p.trial = {-log10(p)}; % change of p value scaling to facilitate reading in topographical plots
        data_p.time = {data.time{1}(cfg.toi)};
        data_p.trialinfo(2:end,:) = [];
        data_p.sampleinfo(2:end,:) = [];
        data_p.label = data.label(cfg.channel);
        % erf
        data_erf = data_p;
        data_erf.trial =  {erf};
    case 'p'
        data_erf = erf;
        data_p = p;
end
