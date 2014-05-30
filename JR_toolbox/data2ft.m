function ft_data = data2ft(cfg,data)
% ft_data= data2ft(cfg,data)
% transform a matrix into a fieldtrip structure
% cfg.fsample
% cfg.time
% cfg.label_start
% nb_chans    = size(data,1);
% nb_samples  = size(data,2);
% nb_trials   = size(data,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg, 'fsample'), cfg.fsample=1;  end
if ~isfield(cfg,'time'),    cfg.time = cumsum(repmat(1/cfg.fsample,1,size(data,2)));end
if ~isfield(cfg,'label_start'),cfg.label_start = 'eeg';                     end

nb_chans    = size(data,1);
nb_samples  = size(data,2);
nb_trials   = size(data,3);

%-- channels
if ~isfield(cfg,'label')
    for chan = 1:nb_chans
        ft_data.label{chan} = [cfg.label_start num2str(chan)];
    end
     ft_data.label =  ft_data.label';
else
    ft_data.label = cfg.label;
end

%-- times
if ~isfield(cfg,'times')
    for trial = 1:nb_trials
        ft_data.time{trial} = cfg.time;
    end
else
    ft_data.time = cfg.times;
end

%-- trials
if ~isfield(cfg,'trial')
    for trial = 1:nb_trials
        ft_data.trial{trial} = squeeze(data(:,:,trial));
    end
else
    ft_data.trial = cfg.trial;
end

%-- fsample
ft_data.fsample = cfg.fsample;

%-- trialinfo
if ~isfield(cfg,'trialinfo')
    ft_data.trialinfo = ones(nb_trials,1);
else
    ft_data.trialinfo = cfg.trialinfo;
end

%-- cfg
if ~isfield(cfg,'cfg')
    ft_data.cfg.description = 'data2ft';
else
    ft_data.cfg = cfg.cfg;
   
end

%-- sample info
if ~isfield(cfg,'sampleinfo')
    ft_data.sampleinfo = [0:nb_trials',1:(nb_trials+1)'] .* nb_samples;
else
    ft_data.sampleinfo = cfg.sampleinfo;
end

%-- hdr
if ~isfield(cfg,'hdr');
    ft_data.hdr.label               = ft_data.label;
    ft_data.hdr.nChans              = nb_chans;
    ft_data.hdr.Fs                  = cfg.fsample;
    ft_data.hdr.nSamples            = max(ft_data.sampleinfo(:,2));
    ft_data.hdr.nSamplesPre         = 0;
    ft_data.hdr.nSamplesnTrials     = 1;
else
    ft_data.hdr = cf.hdr;
end
