function wbicoh(x,frequencies,sampling_frequency)
% x signal
% frequencies: frequencies of interest
% wavelet bicoherence
clear x
for segment = 1:50,x(segment,:) = cos(rand+linspace(-30*pi,30*pi,1000));end
clear x
phase = rand(1,3) * 2*pi;
freqs = [.04 .07];
w = rand;
% classic examples
x = linspace(-100*pi,100*pi,1000);
x(1,:) = sin(2*pi*freqs(1)*x+phase(1)) + sin(2*pi*freqs(2)*x+phase(2)) + ...
    w*(sin(2*pi*(freqs(1)+freqs(2))*x+(phase(1)+phase(2)))) + ...
    (1-w)*sin(2*pi*(freqs(1)+freqs(2))*x+phase(3))+eps;

% segments
x = linspace(-5*pi,5*pi,100);
y = sin(2*pi*freqs(1)*x+phase(1)) + sin(2*pi*freqs(2)*x+phase(2)) + ...
    w*(sin(2*pi*(freqs(1)+freqs(2))*x+(phase(1)+phase(2)))) + ...
    (1-w)*sin(2*pi*(freqs(1)+freqs(2))*x+phase(3))+eps;

signal = []
for seg = 1:50,
    signal = [signal zeros(1,randi(500)) y];
end
x = signal;
plot(x);
% x = ifft(abs(fft(x))); % randomize angles

% nested frequencies
x = cos(linspace(-25*pi,25*pi,10000))*10+ (1+cos(linspace(-25*pi,25*pi,10000)))/2 .* cos(linspace(-500*pi,500*pi,10000));
x = cos(linspace(-30*pi,30*pi,10000)).*sin(linspace(-500*pi,500*pi,10000));
% similar spectrum than nested frequencies
x = cos(linspace(-30*pi,30*pi,10000))+sin(linspace(-500*pi,500*pi,10000));
% x = ifft(abs(fft(x))); % randomize angles
frequencies = 1:size(x,1)/2;
frequencies = 2:1:45;

n_segments = size(x,1);
n_samples = size(x,2);
sampling_frequency = 250;

sigma = 1;

% 0. baseline correct and apply hanning window on each segment
wind = hanning(n_samples);
wind = wind(:);
x = (x - repmat(nanmean(x,2),1,n_samples)) .* repmat(wind',n_segments,1);

% 1. calculate windowed harmonic wavelet transform for each frequency
a = NaN(n_segments,length(frequencies),n_samples);
for segment = 1:n_segments
    fprintf('*');
    for f = 1:length(frequencies)
        a(segment,f,:) = whwt(x(segment,:),frequencies(f),sigma,sampling_frequency);
    end
end

% 2. Calculate wavelet
Pxx = NaN(n_segments,length(frequencies));
for segment = 1:n_segments
    fprintf('*');
    for f = 1:length(frequencies)
        % power spectrum
        Pxx(segment,f) = sum(a(segment,f,:).*conj(a(segment,f,:)));
    end
end
Pxx = Pxx/length(x);

% 3. Calculate wavelet bispectrum (van Milligen et al 1995):
% Bxxx(fj,fk) = integral over TAUS { ax(fj,t)*ax(fk,t)*ax*(fj+fk,t)dt}
Bxxx = NaN(n_segments,length(frequencies),length(frequencies));
% biphase = NaN(n_segments,length(frequencies),length(frequencies),length(x));
for segment = 1:n_segments
    fprintf('*');
    for f1 = 1:length(frequencies)
        for f2 = f1:(length(frequencies)-f1)
            %  bispectrum
            Bxxx(segment,f1,f2) = sum(...
                a(segment,f1,:).*...
                a(segment,f2,:).*...
                conj(a(segment,f1+f2,:)));
            % biphase
%             biphase(segment,f1,f2,:) = angle(a(segment,f1,:)) + angle(a(segment,f2,:)) - angle(a(segment,f1+f2,:));
        end
    end
end
Bxxx = Bxxx/length(x);

% 4. Calculate noramlized squared wavelet bicoherence (WBIC)
norm2 = NaN(n_segments,length(frequencies),length(frequencies)); % squarred norm
for segment = 1:n_segments
    fprintf('*');
    for f1 = 1:length(frequencies)
        for f2 = f1:(length(frequencies)-f1)
            norm2(segment,f1,f2) = ...
                sum(abs(a(segment,f1,:).*a(segment,f2,:)).^2,3) .* ...
                sum(abs(a(segment,f1+f2,:)).^2,3);
        end
    end
end
bxxx2 = abs(Bxxx).^2./(norm2/length(x));


% 5. phase randomization
% classic: B(fj,fk) = E[| X(fj)X(fk)X*(fj+fk) | exp(iR phi_d(fj,+fk))]
% with wavelet: B(fj,fk)  = intergral over TAUS {ax(fj,t) ax(fk,t) ax*(fj+fk,t) exp(iR phi_d(fj,+fk)) dt}
Rand_Bxxx = NaN(100,n_segments,length(frequencies),length(frequencies));
for r = 1:100
    for segment = 1:n_segments
        fprintf('*');
        for f1 = 1:length(frequencies)
            for f2 = f1:(length(frequencies)-f1)
                % power spectrum
                Rand_Bxxx(r,segment,f1,f2) = sum(...
                    squeeze(a(segment,f1,:).*...
                    a(segment,f2,:).*...
                    conj(a(segment,f1+f2,:))).* ...
                    squeeze(exp(1i * rand * biphase(segment,f1,f2,:))));
            end
        end
    end
end
rand_bxxx = sqrt(abs(squeeze(nanmean(Rand_Bxxx))).^2./squeeze((norm2/length(x))));

n_freqs     = 2*length(frequencies);
bic         = zeros(n_freqs,n_freqs);
Pyy         = zeros(n_freqs,1);

% use the hankel mask to fasten the computation through the use of the
% symetrical properties of bicoherence maps
mask = hankel([1:n_freqs],[n_freqs,1:n_freqs-1]);
Yf12 = zeros(n_freqs,n_freqs);

for segment = 1:n_segments
    ys      = x(segment,:)';
    %Yf      = fft(ys)  / n_samples;
    Yf      = Af(segment,:,:);
    Yf      = Yf([frequencies sort(end-frequencies+1)]); % truncate irrelevant frequencies
    CYf     = conj(Yf);
    Pyy     = Pyy + Yf .* CYf;
    Yf12(:) = CYf(mask);
    bic     = bic + (Yf * Yf.') .* Yf12;
end

bic     = bic / n_segments;
Pyy     = Pyy  / n_segments;
mask(:) = Pyy(mask);
bic     = abs(bic).^2 ./ (Pyy * Pyy.' .* mask);
bic     = fftshift(bic) ;
