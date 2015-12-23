function spectra( data, line )
%SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

% Spectrum
figure;
fs = data(line).fs*1e6;
t = (0:double(data(line).SWL)-1)*1/fs;
x = double(data(line).data);%.*exp(-1i*2*pi*5e6*t);
xdft = (1/length(x))*fft(x);
freq = -fs/2:(fs/length(x)):fs/2-(fs/length(x));
plot(freq,abs(fftshift(xdft)));
set(gca,'XTick',[-fs/2:5e6:fs/2-(fs/length(x))]);
set(gca,'XTickLabel',[-fs/2:5e6:fs/2-(fs/length(x))]/1e6);
xlabel('Frequency [MHz]');
ylabel('Magnitude Digital Number');
grid;

end

