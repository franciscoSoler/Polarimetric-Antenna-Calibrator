function [ pwrMTRs ] = prueba( calPaths,inputPower )
% relativePower is the relative power of each RM, 
% input power is the input power and phase in db and degrees, 

%	rel = relativePower( setdiff( 1:size( relativePower,1 ),[1,3] ),: ) - ...
 %		repmat( relativePower( 1,: ),2,1 )
%
a = [-1 5];
b = [-1 0];
c = [-1 0];
    A = [ 2*(1 - a)  	   	 
          2*(a - b)
          2*(b - c)
          2*(1 - c)
          2*(a - c)
          2*(1 - c)
          1 - b
          a - c]

    B = [1 1 0 0 -1 -1 0 0  0  0  0 0  0  0  0  0
         0 0 0 0  0  1 1 0  0 -1 -1 0  0  0  0  0
         0 0 0 0  0  0 0 0  0  0  1 1  0  0 -1 -1
         1 0 1 0  0  0 0 0 -1  0 -1 0  0  0  0  0
         0 0 0 0  0  1 0 1  0  0  0 0  0 -1  0 -1
         1 0 0 1  0  0 0 0  0  0  0 0 -1  0  0 -1
         0 1 0 0  0  0 0 0  0 -1  0 0  0  0  0  0
         0 0 0 0  0  0 1 0  0  0  0 0  0  0 -1  0 ];
	
	%C = ( A'*A )\( A'*B );
%	C = inv(A'*A) * A' * B;
%	att = C*calPaths(:,1) + inputPower;
%	phase = C*calPaths(:,2)
	C = ( A(:,1)'*A(:,1) )\( A(:,1)'*B );
    D = ( A(:,2)'*A(:,2) )\( A(:,2)'*B );

	% las siguientes 2 lineas pueden achicarse a solo 1
	att = C*calPaths(:,1) + inputPower;
	phase = D*calPaths(:,2);
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
