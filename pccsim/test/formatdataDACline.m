function [ out ] = formatDAC(dataNFO, dataDAC ,line)
%FORMATDAC Summary of this function goes here
%   Detailed explanation goes here

datasize = size(dataDAC);
i=1;
%out = nan(datasize(1),1);
    out(1).fs    = dataNFO.fs(line);
    out(1).comp  = dataNFO.comp(line);
    out(1).SWL   = dataNFO.long(line);
    out(1).t     = (0:out(1).SWL-1)/out(1).fs;
    out(1).data  = double(dataDAC(1,1:dataNFO.long(line)));
    
    if(out(i).SWL== 4096)
        out(i).PRF   = 1500;
    elseif(out(i).SWL== 2560)
        out(i).PRF   = 2500;
    end
    
    j=mod(i,0.5*1500+0.5*2500);
    if(j<=0.25*1500) % up prf 1500
        out(i).chirptype = 'up';
        out(i).PRF = 1500;
    elseif(j>0.25*1500 && j<=0.5*1500) % down prf 2500
        out(i).chirptype = 'down';
        out(i).PRF = 1500;
    elseif(j>0.5*1500 && j<=(0.5*1500+0.25*2500))
        out(i).chirptype = 'up';
        out(i).PRF = 2500;
    elseif(j>(0.5*1500+0.25*2500) && j<=(0.5*1500+0.5*2500))
        out(i).chirptype = 'down';
        out(i).PRF = 2500;
    end
end

