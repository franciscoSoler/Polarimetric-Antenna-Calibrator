clear all;
close all;
% clc;
tic
testnb = '43';

% load NFO file
% filename = ['data/l0-acqId', testnb, '-beamId1-polHV_typescience.xml'];
% dataNFO  = readNFO( filename );

% LINES  = 271316; 
LINES  = 251145; 
SAMPLES = 2560;
PRF = 2500;

line_length = SAMPLES * 4;

filename = ['Prefijo_ACQID',testnb,'_M.DACV'];
% filename = ['data/l0-acqId', testnb, '-beamId1-polHV_typescience.bin'];
fid = fopen(filename);
frewind(fid);
signalpower = zeros(LINES, 1);

f0 = 30e6;
fs = 120e6;
l1 = double(SAMPLES/2);
t = (0:SAMPLES-1)'/fs;

for line = 1:LINES
    % load DAC line
    
    if line == 251142
        dataDAC = fread(fid, 1536, 'float32', 'b');
    else
        dataDAC = fread(fid, SAMPLES, 'float32', 'b');
    end
    %{
    dataDAC = fread(fid, SAMPLES, 'float32', 'b');
    %}
    % dataDAC  = readDACline( fid, dataNFO, line );
    
    % subtract mean of a signal (ver posiblemente no se tenga que eliminar)
    dataDAC = detrend(dataDAC')';

    % Spectrum
    % this script obtains the spectrum of the signal, then it keeps with the positive frequencies and then 
    % moves that spectrum to the zero frequency. and returns to time domain. 
    xdft = fftshift((1/SAMPLES)*fft(dataDAC));

    if line == 251142
        Spp = fftshift([zeros(l1, 1); sqrt(2)*xdft(l1+1:1536)]);
        sp = ifft(Spp) .* exp(-1i*2*pi*f0*(0:1536-1)'/fs);
    else
        Spp = fftshift([zeros(l1, 1); sqrt(2)*xdft(l1+1:SAMPLES)]);
        sp = ifft(Spp) .* exp(-1i*2*pi*f0*t);
    end
    %{
    Spp = fftshift([zeros(l1, 1); sqrt(2)*xdft(l1+1:SAMPLES)]);
    sp = ifft(Spp) .* exp(-1i*2*pi*f0*t);
    %}
    %signalpowerpb(line) = var(data.data);
    % Compute and plot spectra
    %spectra(data,line);
    
    % Plot ENV/I/Q signals
    %plotsignals(t,sp,si,sq);
    
    % Signal Power
    signalpower(line) = var(sp);
    % 1233 is one point, it is selected randomly
    signalphase(line) = mod(rad2deg(phase(sp(1233))),360);
    %powbp = 2 * bandpower(s,Fs,[0 Fs/2]);
end

fclose(fid);

toc
time = [1:LINES]'/PRF;
figure;plot(time, 10*log10(signalpower/1e-3));
figure;plot(time, signalphase);

figure;stairs(time, 20*log10(signalpower/1e-3)-max(20*log10(signalpower/1e-3))+9.7);
