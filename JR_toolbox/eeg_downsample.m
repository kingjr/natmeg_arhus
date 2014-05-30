function [dataOut times] = eeg_downsample(data,fmax,srate,times,pnts)
%   EEG_DOWNSAMPLE   Downsample data set based on maximum frequency of interest
%       [EEG] = EEG_DOWNSAMPLE(data,FMAX)
%
%   Uses DCT to minimize artefacts (better than fft)
%
%   Created by Alexandre Gramfort on 2009-04-28.
%   Modified by JR KING 2011-01-04
%   All rights reserved.

if nargin < 2
    fmax = 256; % 256 Hz
    srate = 1024;
    times = [-800:4:700];
    pnts = length(times);
end

resamp = fmax*2/srate;
srate = fmax*2;

disp(['Downsampling dataset to ',num2str(srate),' Hz']);

for e = 1:size(data,1) % electrode
    progressbar(e,size(data,1))
for t=1:size(data,3) % trial
    dataOut(:,:,ii) = dctresample(data(e,:,t)',round(size(data,2)*resamp))';
end
end

times = linspace(times(1),times(end),pnts);

end %  function