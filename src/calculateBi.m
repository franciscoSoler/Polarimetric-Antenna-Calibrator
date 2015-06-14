function B = calculateBi ( xi,yi,qttyAnt )
% inputs: 
%	* xi are the attenuation and phase shift of each compoent, are noted 
% that this vector in this simulation has only 12 rows.
%	* yi corresponds to attenuation and phase shift of the different distances
% between the different RM
%	The format corresponds to an absolute, phase pair

	A = zeros( qttyAnt^2,2*qttyAnt+4 );
	Y = zeros( qttyAnt^2,2 );
	A( :,1:2 ) = 1;
	A( :,size( A,2 )-1:size( A,2 ) ) = 1;
	
	for i = 1:qttyAnt,
		index = ( i-1 ) * qttyAnt;
		inicial = i-1;
		for j = 1:qttyAnt,
			idx = index + j;
			A( idx,2+i ) = 1;
			A( idx,2+qttyAnt+j ) = 1;

			if inicial ~= 0
				Y( idx,: ) = yi( inicial+1,: );
				inicial = inicial -1;
			else
				Y( idx,: ) = yi( j - ( i-1 ),: );
			end
		end
	end
	
	att = A*xi(:,1) + Y(:,1);
	phase = A*xi(:,2) + Y(:,2);
	B = [att phase];
end
