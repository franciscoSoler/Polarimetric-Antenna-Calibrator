clear all;
close all;
clc;

testnb = '43';
LINES  = 2001;
line = 1;

% load NFO file
filename = ['Prefijo_ACQID',testnb,'_M.NFOH'];
dataNFO  = readNFO( filename );

% load DAC file
filename = ['Prefijo_ACQID',testnb,'_M.DACH'];
dataDAC  = readDAC( filename, dataNFO, LINES );

% format DAC file
data = formatdataDAC(dataNFO,dataDAC);

% Compute and plot spectra
spectra(data,line);

% I/Q signals
[sp,si,sq] = pb2bb(data, line);

% Compute and plot spectra
newdata = data;
newdata(line).data = sp;
spectra(newdata,line);



