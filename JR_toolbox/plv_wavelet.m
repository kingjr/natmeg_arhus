function output = plv_wavelet(cfg, data)
% Phase locking value using wavelets (and not hilbert)
% Performs time-frequency decomposition using Fieldtrip and then calculates
% the phase locking value based on phase obtained with complex wavelets
%--------------------------------------------------------------------------
%-- usage : output = plv_wavelet(cfg,data)
%-- input : 
%       + cfg :
%           * foi           : frequencies of interest      
%           * toi           : time of interest in seconds
%           * [method]      : wavelet or fdr               [default : wavelet]
%           * [width]       : wavelet width                [default : 5]
%           * [gwidth]      : length of the used wavelets  [default : 3]
%           * [feedback]    : give feedback on computation [default : true]
%       + data formatted in a fieltrip format ( with fields trial, time, 
%         fsample and label)
%
%-- output : 
%       + output : output: n_chan x n_chan x freq x time (mean across trials)
%
% Note : requires Fieldtrip and signal processing toolboxes
%
%--------------------------------------------------------------------------
 

% parameters
if nargin ~= 2, error('wrong number of arguments'); end
if ~isfield(cfg,{'foi', 'toi'}), error('missing correct cfg fields');end
if ~isfield(data,{'trial','time','fsample', 'label'}), error('missing correct data fields');end
if ~isfield(cfg,'method'),       cfg.method     = 'wavelet';         end
if ~isfield(cfg,'width'),        cfg.width      = 5 ;                end                
if ~isfield(cfg,'gwidth'),       cfg.gwidth     = 3;                 end
if ~isfield(cfg,'feedback'),     cfg.feedback   = true;              end
cfg.output='fourier'; % output of freqanalysis will be complex values
cfg.keeptrials='yes'; % keep trials, so plv can be computed across trials

% time-frequency analysis
freq=ft_freqanalysis(cfg,data); % freq.fourierspctrm: trials x channels x freqs x time

[~, n_chan n_freq n_sample]=size(freq.fourierspctrm);

% create channel pairs matrix
chan_index          = [];%
for ii = 1:(n_chan-1)
    for jj = (ii+1):n_chan
        chan_index(:,end+1) = [ii jj];
    end
end
if rem(n_chan,2)~=0
chan_index          = reshape(chan_index,[2, size(chan_index,2)/n_chan, n_chan]); % reshape for memory issues
else
    chan_index          = reshape(chan_index,[2, size(chan_index,2)/(n_chan-1), (n_chan-1)]); % reshape for memory issues
end

plv_mat=zeros(size(chan_index,2),n_freq,n_sample,size(chan_index,3));
for cut=1:size(chan_index,3)
if cfg.feedback, fprintf('\ncomputing plv for cut %d',cut); end
plv_mat(:,:,:,cut)=squeeze(abs(mean(exp(1i*(...
    angle(freq.fourierspctrm(:,chan_index(1,:,cut),:,:))-angle(freq.fourierspctrm(:,chan_index(2,:,cut),:,:)))))));
end

% reshape output matrix for readability
if cfg.feedback, fprintf('\nreshape\n'); end
plv_mat = reshape(permute(plv_mat,[1 4 2 3]),[], size(plv_mat,2), size(plv_mat,3));

output = NaN(n_chan,n_chan,n_freq,n_sample);
kk = 0;
for ii = 1:(n_chan-1)
    for jj = (ii+1):n_chan
        kk = kk+1;
        output(ii,jj,:,:) = squeeze(plv_mat(kk,:,:));
        output(jj,ii,:,:) = squeeze(plv_mat(kk,:,:));
    end
end