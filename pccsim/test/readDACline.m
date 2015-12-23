function [ data ] = readDACline( fid , dataNFO, LINE)
%READDAC Summary of this function goes here
%   Detailed explanation goes here

data = nan(1,dataNFO.long(LINE));

line = LINE;
status = fseek(fid,dataNFO.off(line)*4,'bof');
if(status==0)
    data(1,1:dataNFO.long(line)) = fread(fid,dataNFO.long(line),'float32','b');
else
    disp('Error reading DAC file');
end

end

