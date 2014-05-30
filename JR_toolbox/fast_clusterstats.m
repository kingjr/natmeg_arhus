function stats = fast_clusterstats(x,y,cfg)

cfg = default(cfg,{'neighbour'; 'ft'; 'tl'; 'neighbour';'stats'});
%% transform data into fieldtrip structure
if size(x,3) == 1, x = repmat(x,[1 1 2]); end
cfg.ft      = default(cfg.ft,{'layout', cfg.layout; 'label', cfg.layout.label});
x           = data2ft(cfg.ft,permute(x,[2 3 1]));
x.dimord    = 'chan_time';

%% apply ft transform
cfg.tl      = default(cfg.tl,{'keeptrials','yes'});
x           = ft_timelockanalysis(cfg.tl,x);

%% calculate neighbourhood
cfg.neighbour = default(cfg.neighbour,{'method','distance';'layout', cfg.layout;'neighbourdist',.1});
neighbours    = ft_prepare_neighbours(cfg.neighbour);
% cfg.neighbour = default(cfg.neighbour,{'method','distance';'layout', cfg.layout;'neighbourdist',.05});
% neighbours    = ft_neighbor_nD(cfg.layout.pnt,cfg.layout.label,cfg.neighbour.neighbourdist);

%% compute stats
cfg.stats   = default(cfg.stats,{...
    'neighbours',       neighbours;...
    'layout',           cfg.layout;...
    'method',           'montecarlo';...        % use the Monte Carlo Method to calculate the significance probability
    'statistic',        'indepsamplesT';...     % use the independent samples F-statistic as a measure to evaluate the effect at the sample level
    'correctm',         'cluster';...
    'clusteralpha',     .1;...                  % alpha level of the sample-specific test statistic that will be used for thresholding
    'clusterstatistic', 'maxsum';...            % test statistic that will be evaluated under the permutation distribution.
    'minnbchan',        2;...                   % minimum number of neighbourhood channels that is required for a selected sample to be included in the clustering algorithm (default=0).
    'tail',             0;...                   % -1, 1 or 0 (default = 0); one-sided or two-sided test
    'clustertail',      0;...
    'alpha',            0.05;...                % alpha level of the permutation test
    'numrandomization', 100;...
    'design',           y;...
    'latency',          x.time;...
    'avgovertime',      'yes';...
    'parameter',        'trial'});
cfg.stats.design(y==min(y))   = 1;
cfg.stats.design(y==max(y))   = 2;

stats                   = ft_timelockstatistics(cfg.stats, x);
%% retrieve significancy
try
    stats.sigchan = zeros(1:length(stats.label),1);
    for c = 1:length(stats.posclusters)
        if stats.posclusters(c).prob < cfg.stats.alpha
            stats.sigchan(stats.posclusterslabelmat==c) = 1;
        end
    end
    for c = 1:length(stats.negclusters)
        if stats.negclusters(c).prob < cfg.stats.alpha
            stats.sigchan(stats.negclusterslabelmat==c) = 1;
        end
    end
    stats.sigchan = find(stats.sigchan==1);
catch e
    stats.sigchan=[];                                      % nothing sig.
end


function x = default(x,fields)
if ~iscell(fields), fields = {fields}; end
if size(fields,2) == 1, fields = cat(2,fields,cell(size(fields,1),1));end
for f = 1:size(fields,1)
    if ~isfield(x,fields{f,1})
        eval(['x.' fields{f,1} '=fields{f,2};']);
    end
end