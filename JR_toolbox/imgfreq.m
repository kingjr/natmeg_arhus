function imgfreq(data,varargin)
if nargin == 1, varargin = {} ; end
for ii = 1:2:length(varargin)
    eval([varargin{ii} '=varargin{ii+1};']);
end
f = hilbert(data);
a = angle(f);
p = abs(f);

if ~exist('scale', 'var'), scale = max(p); end
y = cat(3,(cos(a)+1)/2,(sin(a)+1)/2,p./scale);

y = hsv2rgb(y);
imagesc(y)
