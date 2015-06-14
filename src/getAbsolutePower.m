function pwrMTRs = getAbsolutePower( calPaths,inputPower )
% relativePower is the relative power of each RM, 
% input power is the input power and phase in db and degrees, 

%	rel = relativePower( setdiff( 1:size( relativePower,1 ),[1,3] ),: ) - ...
 %		repmat( relativePower( 1,: ),2,1 )
%
	if size( calPaths,1 ) == 9

		A = [ 1 0 0 1 0 0 1 0 0
			  1 0 0 0 1 0 0 1 0
			  1 0 0 0 0 1 0 0 1
			  0 1 0 1 0 0 0 1 0
			  0 1 0 0 1 0 1 0 0 
			  0 1 0 0 0 1 0 1 0
			  0 0 1 1 0 0 0 0 1
			  0 0 1 0 1 0 0 1 0
			  0 0 1 0 0 1 1 0 0 ]

		pwrMTRs = ( A'*A )\( A'*calPaths );

		return
	%{	
		A = [ 2 -2
			  2  4
			  4  2
			  2  1 ];

		A = [ 2 -2
			  0  2
			  2  0
			  1  0 ];
		

		B = [ 1  1  0 -1 -1  0  0  0  0
			  0  0  0  0  1  1  0 -1 -1
			  1  0  1  0  0  0 -1  0 -1
			  0  1  0  0  0  0  0 -1  0 ];
	%}
	else
		A = [ 0 -2  0  	   	 
			  8  4  2  	   	
			 -8 -2 -4    
			  8  2  2  	  
			  0  2 -2 
			  0  0 -2 
			  4  1  1 
			  0  1 -1 ]; 
	%
	%{	
		A = [ 0 
			  8 
			 -8
			  8
			  0
			  0
			  4
			  0];

	%}
		B = [1 1 0 0 -1 -1 0 0  0  0  0 0  0  0  0  0
			 0 0 0 0  0  1 1 0  0 -1 -1 0  0  0  0  0
			 0 0 0 0  0  0 0 0  0  0  1 1  0  0 -1 -1
			 1 0 1 0  0  0 0 0 -1  0 -1 0  0  0  0  0
			 0 0 0 0  0  1 0 1  0  0  0 0  0 -1  0 -1
			 1 0 0 1  0  0 0 0  0  0  0 0 -1  0  0 -1
			 0 1 0 0  0  0 0 0  0 -1  0 0  0  0  0  0
			 0 0 0 0  0  0 1 0  0  0  0 0  0  0 -1  0 ];
	end
	%C = ( A'*A )\( A'*B );
%	C = inv(A'*A) * A' * B;
%	att = C*calPaths(:,1) + inputPower;
%	phase = C*calPaths(:,2)
	C = ( A'*A )\( A'*B );

	% las siguientes 2 lineas pueden achicarse a solo 1
	att = C*calPaths(:,1) + inputPower;
	phase = C*calPaths(:,2);
%{
	att = [ att(1:2,1)
			-sum(att)
			att(3,1) ];
	phase = [ phase(1:2,1)
			  -sum(phase)
			  phase(3,1) ];
%}
	pwrMTRs = [att phase];
end
