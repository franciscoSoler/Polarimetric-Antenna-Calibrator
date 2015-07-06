function [sp] = pb2bb1( data )
%HAPPY Summary of this function goes here
%   Detailed explanation goes here

fs = 120e6;
fc = 30e6;
bw = 50e6;
tp = 20e-6;
swl = 21.34e-6;

t = 0:1/fs:swl-1/fs;   	 					% Tiempo para cada sample
l1 = length(t)/2;
l2 = length(t);

% Spectrum
% this script obtains the spectrum of the signal, then it keeps with the positive frequencies and then 
% moves that spectrum to the zero frequency. and returns to time domain. 
Sp = fftshift((1/length(data))*fft(data));
Spp = fftshift([zeros(1,l1),sqrt(2)*Sp(l1+1:l2)]);
sp = ifft(Spp) .* exp(-1i*2*pi*fc*t);