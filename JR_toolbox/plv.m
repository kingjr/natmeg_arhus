function [output cfg phase_hilbert_data] = plv(data,cfg)
% calculates phase locking value across trials for each pair of channels
%--------------------------------------------------------------------------
%[output cfg phase_hilbert_data] = plv(data,cfg)
%--------------------------------------------------------------------------
%-- input: 
%     + data: n_chan x n_sample x n_trial
%     + cfg
%         * time              : timing of each sample                     
%         * fsample           : sampling rate                             
%         * filt_freqband     : frequency band of interest
%         * [time_select]     : desired window of interest                [min max]
%         * [cut_nb]          : number of cuts to avoid memory problem    n_chan/2
%         * [filt_order]      : filter order                              [6 6]
%         * [filt_type]       : filter type                               'but'
%         * [filt_dir]        : filter direction                          'twopass'
%         * [void]            : display feedback                          false                            
%         * [filter_method]   : type of filter (bandpass, high+low pass)  'band'
%         * [mean_across]     : time: 2, trials:3                          2
%-- output:
%         + output            : n_chan x n_chan x n_sample (mean across trials)
%                               phase locking value: PLV = abs(output)
%                               phase lag index:     PLI = sign(angle(output))
%         + cfg
%
%-- requires:
%         + fieldtrip 2011 toolbox
%         + signal processing toolbox
%--------------------------------------------------------------------------        
% Imen El Karoui & Jean-Rémi King 2011, All right reserved (c)
%--------------------------------------------------------------------------
%
% update 
%         + 23/04/2011          JR: changed mean by nanmean
%         + 11/10/2011:1659     Jaco: add unwrap to the angle 
%         + 11/10/2011:1638     JR: add mean_across option, and change
%                               default to mean across time
%         +  5/10/2011:1930     JR:  1. output complexe value instead of abs in
%                               order to compute the phase lag index as
%                               well as the PLV.
%                                   2. include hilbert in function in case
%                                   Matlab fixes its hilbert's bug
% todo
%         + remove time select from this function
%         + improve feedback display
%         + use matlab filtering function instead of fieldtrip
%--------------------------------------------------------------------------


%-- parameters
if nargin ~= 2, error('wrong number of arguments'); end
if ~isfield(cfg,{'time', 'fsample', 'filt_freqband'}), error('missing correct cfg fields');end
if ~isfield(cfg,'time_select'),     cfg.time_select     = [min(cfg.time) max(cfg.time)];end
if ~isfield(cfg,'cut_nb'),          cfg.cut_nb          = size(data,1)/2;               end
if ~isfield(cfg,'filt_order'),      cfg.filt_order      = [4 4];                       end
if ~isfield(cfg,'filt_type'),       cfg.filt_type       = 'but';                        end
if ~isfield(cfg,'filt_dir'),        cfg.filt_dir        = 'twopass';                    end
if ~isfield(cfg,'void'),            cfg.void            = false;                        end
if ~isfield(cfg,'filter_method'),   cfg.filter_method   = 'band';                       end
if ~isfield(cfg,'mean_across'),     cfg.mean_across     = 2;                            end % 2:time, 3:tiral
    
%--- select timing
data                = data(:,find(cfg.time>=cfg.time_select(1),1):find(cfg.time>=cfg.time_select(2),1),:);
[n_chan n_sample n_trial] = deal(size(data,1),size(data,2),size(data,3));

%-- channel combination vector
if cfg.void, fprintf('create combination vector'); end
chan_index          = [];%
for ii = 1:(n_chan-1)
    for jj = (ii+1):n_chan
        chan_index(:,end+1) = [ii jj];
    end
end
if size(chan_index,2)/cfg.cut_nb ~= round(size(chan_index,2)/cfg.cut_nb),
    error('cut should divide into integer number of dimensions');
end
chan_index          = reshape(chan_index,[2, size(chan_index,2)/cfg.cut_nb, cfg.cut_nb]); % reshape for memory issues

%-- high & low pass, hilbert angle
if cfg.void, fprintf('\ncomputes filtering, hilbert & angle\n'); end
phase_hilbert_data  = NaN(n_chan, n_sample,n_trial);
parfor trial = 1:n_trial
    if cfg.void, fprintf('+'); end
    switch cfg.filter_method
        %  /!\  hilbert is passed across wrong dimensions!
        case 'high_low'
            phase_hilbert_data(:,:,trial) = angle(hilbert(...
                ft_preproc_lowpassfilter(...
                ft_preproc_highpassfilter(...
                squeeze(data(:,:,trial)),...
                cfg.fsample,cfg.filt_freqband(2),cfg.filt_order(2),cfg.filt_type,cfg.filt_dir),...
                cfg.fsample,cfg.filt_freqband(1),cfg.filt_order(1),cfg.filt_type,cfg.filt_dir)')');
        case 'band'
            phase_hilbert_data(:,:,trial) = angle(hilbert(...
                ft_preproc_bandpassfilter(...
                squeeze(data(:,:,trial)),...
                cfg.fsample,cfg.filt_freqband,cfg.filt_order(1),cfg.filt_type,cfg.filt_dir)')');
    end
end

%keyboard

%-- phase difference
if cfg.void, fprintf('\ncomputes PLV\n'); end
switch cfg.mean_across
    case 3
        plv_mat = NaN(n_sample,size(chan_index,2),cfg.cut_nb);                      % initialize
        output = NaN(n_chan,n_chan,n_sample);
    case 2
        plv_mat = NaN(n_trial,size(chan_index,2),cfg.cut_nb);                      % initialize
        output = NaN(n_chan,n_chan,n_trial);
    otherwise
        error('wrong dim in cfg.mean_cross');
end
parfor cut = 1:cfg.cut_nb                                                   % cut vector to avoid massively big matrices
    if cfg.void, fprintf('+'); end
    phase_diff = unwrap(phase_hilbert_data(chan_index(1,:,cut),:,:))-unwrap(phase_hilbert_data(chan_index(2,:,cut),:,:));
    plv_mat(:,:,cut) = squeeze(nanmean(exp(1i*phase_diff),cfg.mean_across))';  % 2 = mean PLV across time; 3 = mean across trials;
end

%-- reshape for matrix readibility
if cfg.void, fprintf('\nreshape'); end
plv_mat = reshape(plv_mat, size(plv_mat,1), []);


kk = 0;
for ii = 1:(n_chan-1)
    for jj = (ii+1):n_chan
        kk = kk+1;
        output(ii,jj,:) = plv_mat(:,kk);
    end
end

function x = hilbert(xr,n)
% duplicate matlab hilbert function in case they fix their bug
% ------------------------------------------------------------------------
%HILBERT  Discrete-time analytic signal via Hilbert transform.
%   X = HILBERT(Xr) computes the so-called discrete-time analytic signal
%   X = Xr + i*Xi such that Xi is the Hilbert transform of real vector Xr.
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%   If Xr is a matrix, then HILBERT operates along the columns of Xr.
%
%   HILBERT(Xr,N) computes the N-point Hilbert transform.  Xr is padded with 
%   zeros if it has less than N points, and truncated if it has more.  
%
%   For a discrete-time analytic signal X, the last half of fft(X) is zero, 
%   and the first (DC) and center (Nyquist) elements of fft(X) are purely real.
%
%   EXAMPLE:
%          Xr = [1 2 3 4];
%          X = hilbert(Xr)
%          % produces X=[1+1i 2-1i 3-1i 4+1i] such that Xi=imag(X)=[1 -1 -1 1] is the
%          % Hilbert transform of Xr, and Xr=real(X)=[1 2 3 4].  Note that the last half
%          % of fft(X)=[10 -4+4i -2 0] is zero (in this example, the last half is just
%          % the last element).  Also note that the DC and Nyquist elements of fft(X)
%          % (10 and -2) are purely real.
%
%   See also FFT, IFFT.

%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.10.4.4 $  $Date: 2008/09/13 07:14:19 $

%   References:
%     [1] Alan V. Oppenheim and Ronald W. Schafer, Discrete-Time
%     Signal Processing, 2nd ed., Prentice-Hall, Upper Saddle River, 
%     New Jersey, 1998.
%
%     [2] S. Lawrence Marple, Jr., Computing the discrete-time analytic 
%     signal via FFT, IEEE Transactions on Signal Processing, Vol. 47, 
%     No. 9, September 1999, pp.2600--2603.

if nargin<2, n=[]; end

if ~isreal(xr)
  warning(generatemsgid('Ignore'),'HILBERT ignores imaginary part of input.')
  xr = real(xr);
end
% Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr);
if isempty(n)
  n = size(xr,1);
end
x = fft(xr,n,1); % n-point FFT over columns.
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.
if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
elseif n>0
  % odd and nonempty
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end
x = ifft(x.*h(:,ones(1,size(x,2))));

% Convert back to the original shape.
x = shiftdim(x,-nshifts);