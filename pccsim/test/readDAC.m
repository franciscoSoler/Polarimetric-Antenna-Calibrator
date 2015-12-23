function [ data ] = readDAC( filename , dataNFO, LINES)
%READDAC Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename);
frewind(fid);

data = nan(LINES,dataNFO.long(1));

for line = 1:LINES
    status = fseek(fid,dataNFO.off(line)*4,'bof');
    if(status==0)
        data(line,1:dataNFO.long(line)) = fread(fid,dataNFO.long(line),'float32','b');
    else
        disp('Error reading DAC file');
    end
end
fclose(fid);

end

