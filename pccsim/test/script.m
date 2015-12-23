%
clear all
close all
% clc
tic
testnb = '43';

% load NFO file
filename = ['Prefijo_ACQID',testnb,'_M.NFOV'];
dataNFO  = readNFO( filename );

LINES  = dataNFO.LINESNUMBER;

filename = ['Prefijo_ACQID',testnb,'_M.DACV'];
fid = fopen(filename);
frewind(fid);

for line=1:LINES
    % load DAC line
    
    dataDAC = readDACline( fid, dataNFO, line );
    
    % subtract mean of a signal (ver posiblemente no se tenga que eliminar)
    dataDAC = detrend(dataDAC')';
    
    % format DAC file
    data = formatdataDACline(dataNFO,dataDAC,line);
    
    %signalpowerpb(line) = var(data.data);
    % Compute and plot spectra
    %spectra(data,line);
    
    % I/Q signals
    [t,sp,si,sq] = pb2bb1(data);
    
    % Plot ENV/I/Q signals
    %plotsignals(t,sp,si,sq);
    
    % Signal Power
    signalpower(line) = var(sp);
    signalphase(line) = mod(rad2deg(phase(sp(1233))),360);
    %powbp = 2 * bandpower(s,Fs,[0 Fs/2]);
end


figure;plot((1:length(signalpower))/2500, 10*log10(signalpower/1e-3));
figure;plot((1:length(signalpower))/2500,signalphase);

figure;stairs((1:length(signalpower))/2500,20*log10(signalpower/1e-3)-max(20*log10(signalpower/1e-3))+9.7);

% select steps
data = 10*log10(signalpower/1e-3);

delta = 0.5;
steps = 32;
starts = zeros(steps, 1);
ends = zeros(steps, 1);
starts(1) = 121655;
ends(steps) = 250039;
j = 1;
notErased = true;
for i = 121655:250034 % tengo 4 para irme para adelante 
    if abs(data(i) - data(i+1)) > delta || (i - starts(j)) / 2500 > 2
        if j == 4 && notErased
            starts(j) = i+1;
            ends(j - 1) = i;
            notErased = false;
        else
            starts(j+1) = i+1;
            ends(j) = i;
            j = j + 1;
        end
    end
end

% acomodate order pulses
starts = [starts(6:32); starts(1:5)];
ends = [ends(6:32); ends(1:5)];

% select pulses
amountPulses = 10;
startPulses = floor((ends - starts - amountPulses)/2) + starts;
in = zeros(steps * amountPulses, 2560);

fid = fopen(filename);
frewind(fid);
SAMPLES = 2560;

for i = 0:steps-1
    for j = 1: amountPulses
        % in(i*amountPulses + j) = chirpBp(startPulses(i+1) + j - 1, :);
        fseek(fid, (startPulses(i+1) + j - 1)* SAMPLES*4,'bof');
        in(i*amountPulses + j, :) = fread(fid, SAMPLES, 'float32', 'b');
    end
end
fclose(fid);
%
% generate PCCs
fs = 120e6;
fc = 0;
bw = 50e6;
tp = 20e-6;
swl = 21.34e-6;
%chirp = createChirp( fs,fc,bw,tp,swl );
%chirp = load('ce_loop_cr');
chirp = load('AntennaChirpReplica');


walsh_data = load('data');
walsh = walsh_data.walsh_phi_m;

gain = zeros(20, 10);
phase = zeros(20, 10);
for i = 1:amountPulses
    % [gain(:, i) phase(:, i)] = pccdecoder(in(i:amountPulses:length(in(:,1)), :), chirp, walsh);
    % [gain(:, i) phase(:, i)] = pccdecoder(in(i:amountPulses:length(in(:,1)), :), pb2bb(chirp.chirpreplica(10,:)), walsh);
    [gain(:, i) phase(:, i)] = pccdecoder(in(i:amountPulses:length(in(:,1)), :), pb2bb(chirp.dataDAC), walsh);
end
toc


for i = 0:steps-1
    for j = 1: amountPulses
        toPrint(i*amountPulses + j) = data(startPulses(i+1) + j - 1);
    end
end
figure; spectrogram(real(createChirp( fs,30e6,bw,tp,swl )), 128,120,1024,120e6);
figure; spectrogram(chirp.chirpreplica(10,:), 128,120,1024,120e6);