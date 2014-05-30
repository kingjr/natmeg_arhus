function ftdata = raw2fieldrip(data,cfg)
% takes raw data: channels * samples * trial
% % output fieldtrip structure
% cfg.Fs                % sampling frequency
% cfg.nSamplesPre       % pretrigger sample
% cfg.label             % labels
% cfg.orig           	
% cfg.time              %[min(time), max(time)]


if nargin == 1
    cfg.tmp = [];
end

%-- make labels
labels = {};
for ii = 1:size(data,1);
    if ii < 10
        nb = ['00' num2str(ii)];
    elseif ii < 100
        nb = ['0' num2str(ii)];
    else
        nb = [num2str(ii)];
    end
    labels{length(labels) + 1,1} = ['EEG' nb];
end


ftdata.hdr.Fs             = getoption(cfg,'Fs',250);        % sampling frequency
ftdata.hdr.nChans         = size(data,1);                   % number of electrodes
ftdata.hdr.nSamples       = size(data,2);                   % number of samples per trial
ftdata.hdr.nSamplesPre    = getoption(cfg,'nSamplesPre',0); % pretrigger sample
ftdata.hdr.nTrials        = size(data,3);
ftdata.hdr.label          = getoption(cfg,'chanLabels',labels);
ftdata.hdr.orig           = getoption(cfg,'orig',{});

ftdata.label              = ftdata.hdr.label;
ftdata.time               = getoption(cfg,'time', [0 size(data,2)]);
ftdata.time               = repmat({ftdata.time}, 1, size(data,3));
ftdata.fsample            = ftdata.hdr.Fs;
ftdata.cfg.trl            = NaN .* ones(ftdata.hdr.nTrials,3);

for trial = 1:size(data,3)
ftdata.trial{trial}       = squeeze(data(:,:,trial));
end
return