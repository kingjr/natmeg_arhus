function rgb = angle2rgb(data,mM)
% rgb = angle2rgb(data,mM)
% input complex (2D) data

% scale
if nargin==1, mM = 1./[min(data(:)) max(data(:))]; end

% handle n-dimensional matrices
sz = size(data);
data = reshape(data,[],1);

%
[th r] = cart2pol(real(data),imag(data));
r(r(:)<mM(1)) = mM(1); % min plot value
r(r(:)>mM(2)) = mM(2); % max plot value
r = (r-min(r(:)))./(max(r(:))-min(r(:))); % value between 0 and 1
[data(:,1) data(:,2)] = pol2cart(th,r);
data = data(:,1)+1i*data(:,2);

% transform
hue         = (pi+angle(data))/(2*pi);
saturation  = abs(data);
rgb         = hsv2rgb(cat(2,hue,saturation,ones(length(data),1)));

rgb = reshape(rgb,[sz 3]);