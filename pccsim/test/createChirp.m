function [chirp] = createChirp( fs,fc,bw,tp,swl,varargin )
% This function creates a chirp, depending on the inputs it can contain 
% errors
% Inputs:
%	fs: frequency sampling
%	fc: central frequency
%	bw: bandwidth 
%	tp: pulse duration
%	swl: sampling windows duration
%	varargin: contains the data referent to a chirp with errors
%		attenuation: the word attenuation, it is assumed that the next number
%					contains the attenuation 	
%		att: the attenuation in db
%		phaseShift: the word phaseShift, it is assumed that the next number
%					contains the phaseShift in deg
	phShift = 0;
	att = 1;
	if length( varargin )
		if rem( length( varargin ),2 ) ~= 0
			error( 'the quantity of inputs must be pair' );
		end

		for i = 1:2:length( varargin )
			if strcmp( varargin{ i },'phaseShift' )
				phShift = varargin{ i+1 } * pi/180; %[rad] 
			elseif strcmp( varargin{ i },'attenuation' )
				att = db2mag( -varargin{ i+1 } ); % [V]
			else
				error( '%s is not accepted as input',varargin{i} ); 
			end
		end
	end

	kr = bw/tp;                     			% Chirp rate en rango
	t = 0:1/fs:swl-1/fs;   	 					% Tiempo para cada sample
	f0 = fc-bw/2;                   			% Frecuencia inicial
	chirp = att * exp( 1i*phShift ) * exp( 1i*2*pi*( f0*t + kr/2*t.^2 ) ); % chirp ideal, con envolvente unitaria
end
