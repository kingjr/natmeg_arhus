function trl = ft_jr_art_trialfun(cfg)
%% trl = ft_jr_art_trialfun(cfg)
% function that automatically detect the ECG beats, and gives back the fieldtrip trial structure
% - inputs
%   - cfg.
% - output
%   - trl structure (matrix n trials x 3 colums)
% -------------------------------------------------------------------------
% requires 
%   - fieldtrip 2011
%   - (MNE toolbox)
% -------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) 2011 Jean-Rémi KING
% jeanremi.king+matlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global threshold data_artchan
%-- trialfun function used in ft_definetrial to realign the ecg to the
% heart beat
%-- requires cfg.dataset
if ~isfield(cfg, 'dataset'),    error('needs cfg.dataset! (file path)');        end
if ~isfield(cfg, 'artchan'),    error('needs cfg.artchan! (artefact channel)'); end
if ~isfield(cfg, 'triggerchan'),cfg.triggerchan = false; end
if ~isfield(cfg, 'art_prestim'),cfg.prestim     = .100;                         end % select time before heart beat
%if ~isfield(cfg, 'art_poststim'),cfg.poststim    = .550;                         end % select time after heart beat
if ~isfield(cfg, 'threshold'),  cfg.threshold   = 3.5;                          end % in STD
if ~isfield(cfg, 'continuous'), cfg.continuous  = 'yes';                        end % continuous data
if ~isfield(cfg, 'hpfilter'),   cfg.hpfilter    = 'yes';                        end % remove shifts
if ~isfield(cfg, 'hpfreq'),     cfg.hpfreq      = 2;                            end % high pass filter freq
if ~isfield(cfg, 'maxartfreq'), cfg.maxartfreq  = .150;                         end % max s between heart beat 
if ~isfield(cfg, 'minnb'),      cfg.minnb       = 5;                            end % minimum number of artefact
if ~isfield(cfg, 'triggerchan'),cfg.triggerchan = false;                        end
if ~isfield(cfg, 'derivative'), cfg.derivative  = false;                        end 
if ~isfield(cfg, 'lag'),        cfg.lag         = .02;                          end % n lag derivative in s
if ~isfield(cfg, 'refractory'), cfg.refractory  = .5;                           end % refractory period
if ~isfield(cfg,'artmethod'),   cfg.artmethod = 'abs_mean'; end

%% Define generic functions
%------ automatic printing
if cfg.void
    print = @(x) disp([repmat('-',1,20) ' ' x ' ' repmat('-',1,50-length(x))]);
else
    print = @(x) false;
end

%% Main
%-- select ecg chan only
cfg_art                     = cfg;
cfg_art.channel             = cfg.artchan(:);
%-- load data
print('Reads continuous data');
cfg.detrend              = 'yes';
evalc('data_ecg             = ft_preprocessing(cfg_art);');
while 1
    fsample=data_ecg.fsample;
    %-- build dipole
    if length(data_ecg.label) == 1
        data_artchan = data_ecg.trial{1};
    elseif mod(length(data_ecg.label) == 1,2)==0
        switch cfg.artmethod
            case 'subtract',data_artchan = nanmean(data_ecg.trial{1}(1:2:end,:)-data_ecg.trial{1}(2:2:end,:),1);
            case 'mean',    data_artchan = nanmean(data_ecg.trial{1},1);
            case 'abs_mean',data_artchan = nanmean(abs(data_ecg.trial{1}));
        end
    else
        error('the number of artefacted channels should either be 1 or a multiple of 2');
    end
    %-- whiten
    data_artchan = data_artchan - median(data_artchan);
    if cfg.derivative
        %-- compute n-lag derivative
        lag          = ceil(cfg.lag*fsample);
        data_artchan = smooth(data_artchan,lag)';
        data_artchan = abs(data_artchan(:,lag:end)-data_artchan(:,1:(end-lag+1)));
        data_artchan(end+1) = 0;
    end
    %-- find pulse
    print('Finds threshold crossing');
    data_ecg.pulse              = [];
    threshold                   = (median(data_artchan)+cfg_art.threshold*mad(data_artchan));
    data_ecg.pulse(1,:)         = find(abs(data_artchan) > threshold);
    %-- remove consecutive points: take only points that are further away
    
    %-- remove consecutive trials
    bad_trials          = (data_ecg.pulse(2:end)-data_ecg.pulse(1:(end-1)))<(cfg.refractory*fsample);
    data_ecg.pulse(bad_trials)     = [];
    
    %-- time periods excluded by the user
    if isfield(cfg,'rmsamples'),    
        data_ecg.pulse          = data_ecg.pulse(~ismember(data_ecg.pulse,cfg.rmsamples));
    end

    
    if cfg_art.threshold < 1
        warning(['!PROBLEM! NO ARTEFACT FOUND']);
        threshold = NaN;
        trl = NaN;
        return
    elseif length(data_ecg.pulse) >= cfg.minnb 
        break
    else
        cfg_art.threshold = cfg_art.threshold * 2/3;
        warning(['! NO ARTEFACT FOUND: => diminish threshold to ' num2str(cfg_art.threshold)]);
    end
end

%-- calculate post stim according to median time across beats
if ~isfield(cfg, 'poststim') || strcmp(cfg.poststim, 'auto')
    cfg.poststim = median(data_ecg.pulse(2:end)-data_ecg.pulse(1:(end-1))) / data_ecg.fsample - cfg.prestim ;
end
%-- final triggers
trl = [...
    ceil(data_ecg.pulse-(cfg.prestim * data_ecg.fsample)); ...
    floor(data_ecg.pulse+(cfg.poststim * data_ecg.fsample)); ...
    floor(repmat(-cfg.prestim .* data_ecg.fsample,1,length(data_ecg.pulse)))]'; 
%-- find trigger after first trigger
if cfg.triggerchan
    print('Reads artchan continuous data');
    cfg_trigger         = cfg_art;
    cfg_trigger.channel = cfg.triggerchan;
    cfg_trigger.hpfilter = 'no';
    evalc('data_trigger                    = ft_preprocessing(cfg_trigger);');
    print('Find first and last triggers');
    triggers = find(round(data_trigger.trial{1})>0);
    print('Remove trials outside the latter');
    trl = trl(intersect(...
        find(trl(:,1)>(triggers(1)-cfg.prestim*data_ecg.fsample)),...
        find(trl(:,2)<(triggers(end)+cfg.poststim*data_ecg.fsample))),:);
end

%-- remove artefact at beginning and end
print('Remove trials outside recording');
trl = trl(intersect(find(trl(:,1)>(cfg.prestim*data_ecg.fsample)), find(trl(:,2)<(data_ecg.sampleinfo(2)-cfg.poststim*data_ecg.fsample))),:);
clf;plot(data_artchan);hold on; scatter(trl(:,1),threshold*ones(size(trl,1),1),'r')
end
