function  data  = readNFO( filename )
%READNFO Summary of this function goes here
%   Detailed explanation goes here

%# type and size in byte of the record fields
recordType = {'uint8' 'uint8' 'uint16' 'uint64'};
recordLen = [1 1 2 8];
R = cell(1,numel(recordType));

%# read column-by-column
fid = fopen(filename,'rb');
for i=1:numel(recordType)
    %# seek to the first field of the first record
    fseek(fid, sum(recordLen(1:i-1)), 'bof');

    %# % read column with specified format, skipping required number of bytes
    R{i} = fread(fid, Inf, ['*' recordType{i}], sum(recordLen)-recordLen(i),'b');
end
fclose(fid); 
clearvars -except R testcs

data.comp  = double(R{1});
data.fs    = double(R{2});
data.long  = double(R{3});
data.off   = double(R{4});
data.LINESNUMBER = length(data.comp);
data.SWLPULSES = unique(data.long);
data.LINESTOPSAR = unique(diff(find(diff(data.long))));

end

