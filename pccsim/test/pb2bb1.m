function [t, sp,si,sq ] = pb2bb( data )
%HAPPY Summary of this function goes here
%   Detailed explanation goes here

f0 = 30e6;
l1 = double(data(1).SWL)/2;
l2 = double(data(1).SWL);

% Spectrum
% this script obtains the spectrum of the signal, then it keeps with the positive frequencies and then 
% moves that spectrum to the zero frequency. and returns to time domain. 
fs = data(1).fs*1e6;
t = (0:double(data(1).SWL)-1)*1/fs;
x = data(1).data;
xdft = (1/length(x))*fft(x);

Sp = fftshift(xdft);
Spp = fftshift([zeros(1,l1),sqrt(2)*Sp(l1+1:l2)]);
spp = ifft(Spp);
sppc = spp .* exp(-1i*2*pi*f0*t);
sp = sppc;
si = real(sp);
sq = imag(sp);


% figure;
% title(['Line ',num2str(1)]);
% subplot(2,1,1);plot(t,si);
% xlabel('Time [\mu s]')
% ylabel('Level');
% grid;
% axis([0 60e-6 -0.5 0.5]);
% set(gca,'XTick',0:5e-6:60e-6);
% set(gca,'YTick',-0.5:0.1:0.5);
% set(gca,'XTickLabel',0:5:60);
% subplot(2,1,2);plot(t,sq);
% xlabel('Time [\mu s]')
% ylabel('Level');
% grid;
% axis([0 60e-6 -0.5 0.5]);
% set(gca,'XTick',0:5e-6:60e-6);
% set(gca,'YTick',-0.5:0.1:0.5);
% set(gca,'XTickLabel',0:5:60);
end

