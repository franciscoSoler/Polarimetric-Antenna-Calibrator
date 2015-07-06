function [ walkm ] = createPhWalshMtx( varargin )
% Inputs
%	k: orden de la secuencia de walsh
%	phShift: en rad
	if not( length( varargin ) )
		error( 'there must be at least two parameters' );
	end
	
	if rem( length( varargin ),2 ) ~= 0
		error( 'the quantity of inputs must be pair' );
	end

	phErr = 0;
	k = 0;
	phShift = 0;
	qttyRows = 1;
	walshEntered = false;
	for i = 1:2:length( varargin )
		if strcmp( varargin{ i },'phErr' )
			phErr = varargin{ i+1 }; % [rad]
		elseif strcmp( varargin{ i },'order' )
			k = varargin{ i+1 }; 
		elseif strcmp( varargin{ i },'phShift' )
			phShift = varargin{ i+1 }; 
		elseif strcmp( varargin{ i },'walshPhMtx' )
			walkm = varargin{ i+1 }; 
			walshEntered = true;
		elseif strcmp( varargin{ i },'duplicateRows' )
			qttyRows = varargin{ i+1 };
		else	
			error( '%s is not accepted as input',varargin{i} ); 
		end
	end

	if not( walshEntered )
		walkm = walsh_mtx(k) * phShift;
	end
	walkm = walkm + randn(size(walkm)) * phErr;
	
	
	walkm = walkm(ceil((1:qttyRows*size(walkm,1))/qttyRows),:);
	
end
