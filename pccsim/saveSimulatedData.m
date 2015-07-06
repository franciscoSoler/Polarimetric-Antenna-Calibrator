function saveSimulatedData( file,mode,data,lenRow,fillValue )
% this function save the data column by column
% inputs:
%	file: is the fileName
%	mode: corresponds to w or a (write or append)
%	data: is the data to save
%	lenRow: corresponds to the line length
%	fillValue : is the value of the filled data
	% Row prefix bytes (2)
	lenPrefix = 2;
	[dataSize, qttyRows] = size( data );	
	large = typecast( uint16( dataSize ),'uint8' )';
	prefix = repmat( large( end:-1:1 ),1,qttyRows );

	%filledData
	filledData = ones( lenRow-lenPrefix-dataSize,qttyRows )*fillValue;
	% Row prefix + data + filledData
	data2save = [ prefix;data;filledData ];

	fid = fopen( file,mode );
	fwrite( fid,round( data2save ),'int8','l' );
	fclose( fid );

end
