function [psi, stdpsi, psisum, stdpsisum]=data2psi2(data,cfg)
% calculates phase slope index 
%
% [psi, stdpsi, psisum, stdpsisum]=data2psi(data [,cfg]);
%
% Input:
%       data:       chan x sample x trial
%       cfg.freqbins:  KxQ matrix. Each row contains the frequencies (in bins), 
%                   over which  PSI is calculated. (freqbins includes the 
%                   last frequency (f+delta f), i.e. the band F in the paper 
%                   is given for the k.th row as F=freqbins(k,1:end-1).
%                   By setting freqbins=[] PSI is calculated across all 
%                   frequencies (wide band).
%       cfg.nepochjack: number of iteration for estimation of std (default:
%                   number of trials)
%       cfg.method:    estimation method (default: jacknife)
%
% Output:
%       psi:       non-normalized PSI values. For M channels PSI is either 
%                   an MxM matrix (if freqbins has one or zero rows)
%                   or an MxMxK tensor if freqbins has K rows (with K>1).
%                   psi(i,j) is the (non-normalized) flow from channel i to
%                   channel j, (e.g., channel i is the sender if psi(i,j) is
%                   positive.)
%       stdpsi:    estimated standard deviation for PSI.
%                   PSI in the paper is given by psi./(stdpsi+eps) (eps is 
%                   included to avoid 0/0 for the diagonal elements)
%       psisum:    sum(psi,2) is the net flux for each channel.
%       stdpsisum: is the estimated standard deviation of psisum. 
%                   (stdpsisum cannot be calculated from psi and stdpsi 
%                   therefore the extra output)
%
% Can work with parallel toolbox after launching matlabpool
%
% Requires 
%               multiprod.m & multitransp.m available at:
%               http://www.mathworks.com/matlabcentral/fileexchange/8773/
%
% adapted from Nolte G, Ziehe A, Nikulin VV, Schlögl A, Krämer N, Brismar
% T, Müller KR.
%    Robustly estimating the flow direction of information in complex physical systems.
%    Physical Review Letters 2008.
%    (for further information:     http://ml.cs.tu-berlin.de/causality/ )
%
% (c) JeanRémi KING 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- set default parameters
[nchan ndat nep]                                = size(data);               % dimensionalities
if nargin == 1,                 cfg             = [];           end         % default parameters
if ~isfield(cfg, 'freqbins'),   cfg.freqbins    = [];           end         % freqs to be analyzed (in bins?)
if ~isfield(cfg, 'method'),     cfg.method      = 'jackknife';  end         % std estimation method
if ~isfield(cfg, 'nepochjack'), cfg.nepochjack  = nep;          end         % number of iteration for std estimation (default: number ofd trials)
if isempty(cfg.freqbins)
    maxfreqbin  = floor(ndat/2)+1;
    cfg.freqbins= 1:maxfreqbin;
else
    maxfreqbin  = max(max(cfg.freqbins));
end
freqnb          = size(cfg.freqbins,1);

%-- main script
cs              = data2cs_event(data,maxfreqbin);                           % main equation
psall           = zeros(nchan,nchan,freqnb);                                % initialize
pssumall        = zeros(nchan,freqnb);
parfor freq=1:freqnb;                                                       % parallel computation possible
    psall(:,:,freq)=cs2ps(cs(:,:,cfg.freqbins(freq,:)));
    pssumall(:,freq)=sum(psall(:,:,freq),2);
end
psi             = squeeze(psall);
psisum          = squeeze(pssumall);

%-- std estimation
csall           = cs;
psloc           = zeros(nchan,nchan,cfg.nepochjack,freqnb);
switch cfg.method
    case 'jackknife'
        for epoch=1:cfg.nepochjack
            csloc       = data2cs_event(data(:,:,epoch),maxfreqbin);
            cs          = (cfg.nepochjack*csall-csloc)/(cfg.nepochjack+1);
            parfor freq=1:freqnb;                                           % parallel computation possible
                psloc(:,:,epoch,freq)=cs2ps(cs(:,:,cfg.freqbins(freq,:)));
            end
        end
        stdpsi          = squeeze(std(psloc,0,3))*sqrt(cfg.nepochjack);
        stdpsisum       = squeeze(std(squeeze(sum(psloc,2)),0,2))*sqrt(cfg.nepochjack);
    case 'none'
        stdpsi          = 0;
        stdpsisum       = 0;
end
return

function ps=cs2ps(cs)
% /!\ to be optimized 
df              = 1;  % /!\ I don't know what this is ...
pp              = cs; % initialize
for f=1:size(cs,3)
    pp(:,:,f)   = diag(cs(:,:,f))*diag(cs(:,:,f))';
end
pp              = cs./ sqrt(pp);
ps              = sum(imag(conj(pp(:,:,1:end-df)).*pp(:,:,1+df:end)),3);
return

function cs=data2cs_event(data,maxfreqbin)
% cs=data2cs_event(data,maxfreqbin)
%
% calculates cross-spectra from data for event-related measurement
% input:
%               data:       chan x sample x trial
%               maxfreqbin: max frequency in bins
% output:
%               cs:         nchan by chan by maxfreqbin by nseg tensor 
%                           cs(:,:,f,i)  contains the cross-spectrum at 
%                           frequency f and segment i
%
% main equation from Nolte et al 2008, Physical Review letters
% NB: sum(conj) = conj(sum)
% cs = conj(sum((prod a a' for each freq and each epoch), over epochs (dim 4)))
% (c) Jean-Rémi King 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data                = permute(data,[2 3 1]);                                % FIXME: change dims form start 
[ndat nep nchan]    = size(data);
maxfreqbin          = min([maxfreqbin,floor(ndat/2)+1]);
%-- hanning window
mywindow            = repmat(hanning(ndat),1,nchan);
%-- fast fourier
datalocfft          = fft(data.*permute(repmat(mywindow,[1 1 size(data,2)]), [1 3 2]));
%-- organize matrix for fast multiplication (get to be mulitplied dim in first two dims)
datalocfft2         = reshape(permute(datalocfft(1:maxfreqbin,:,:),[3 1 2]), [size(datalocfft,3) 1 maxfreqbin size(datalocfft,2)]);
%-- main equation:
cs                  = conj(sum(multiprod(conj(datalocfft2),multitransp(datalocfft2)),4))/nep;
%av                  = squeez(sum(conj(datalocfft(1:maxfreqbin,:,:)),2))/np;% to be check
return