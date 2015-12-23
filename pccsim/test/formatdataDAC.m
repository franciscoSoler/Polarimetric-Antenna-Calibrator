function [ out ] = formatDAC(dataNFO, dataDAC )
%FORMATDAC Summary of this function goes here
%   Detailed explanation goes here

datasize = size(dataDAC);

%out = nan(datasize(1),1);
for i=1:datasize(1)
    out(i).fs    = dataNFO.fs(i);
    out(i).comp  = dataNFO.comp(i);
    out(i).SWL   = dataNFO.long(i);
    out(i).t     = (0:out(i).SWL-1)/out(i).fs;
    out(i).data  = double(dataDAC(i,1:dataNFO.long(i)));
    
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
end

