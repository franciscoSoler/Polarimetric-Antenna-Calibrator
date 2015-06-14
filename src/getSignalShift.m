function [bTransCal bRecCal] = getSignalShift( calPaths )
% calPaths corresponds to attenuation and phase shift of each path into the 
% antenna

	% qtty of paths 
	b = size( calPaths,1 );
	% qtty of antennas
	n = sqrt(b);
	b = ceil(n^3/12 + 3*n^2/8 - 7*n/12);

%	if mod( n,2 ) ~= 0
%		n = 2/3 * (1 + sqrt(1 + 3*b - 3/4));
%	end
	% the second index, representing the paths equally spaced
	iSeg =	b - ceil(n*(n-2)/4);
%	iSeg1 = n*(n-1)/2

	A = zeros( b,n-1 );
	BTrans = zeros( b,n^2 );
	BRec = zeros( b,n^2 );
%{
	% aca elimino la 3º columna
	A = [ 2 -2  0 	   	 
		  2  4  2 	
		 -2 -2 -4   
		  4  2  2 
		  0  2 -2
		  2  0 -2
		  2  1  1
		  0  1 -1 ];
%}	
	row = 0;
	countEq = 0;
	
	for diff = 1:n-1,
		for i = 1:n-diff,
			%row = i + count;
			% first paths 
			for j = 0:(diff-1)/2,
				row = row + 1;
				if i == 1
					A( row,i-1+diff ) = -2;
					A( row,: ) = A( row,: ) - 2;
				else
					A( row,i-1 ) = 2;
					A( row,i-1+diff ) = -2;
				end
				BTrans( row,1+(i-1)*(n+1)+j ) = 1;
				BTrans( row,1+(i-1)*(n+1)+diff-j ) = 1;
				BTrans( row,1+(i+diff-1)*(n+1)-diff+j ) = -1;
				BTrans( row,1+(i+diff-1)*(n+1)-j ) = -1;
				
				BRec( row,1+(i-1)*(n+1)+j*n ) = 1;
				BRec( row,1+(i-1)*(n+1)+(diff-j)*n ) = 1;
				BRec( row,1+(i+diff-1)*(n+1)-(diff-j)*n ) = -1;
				BRec( row,1+(i+diff-1)*(n+1)-j*n ) = -1;
			end
			% paths equally spaced
			if mod(diff,2) == 0
				idx = iSeg + countEq + i;
				if i == 1
					A( idx,i-1+diff ) = -1;
					A( idx,: ) = A( iSeg+countEq+i,: ) - 1;
				else
					A( idx,i-1 ) = 1;
					A( idx,i-1+diff ) = -1;
				end
				BTrans( idx,1+diff/2+(n+1)*(i-1) ) = 1;
				BTrans( idx,1-diff/2+(n+1)*(i+diff-1) ) = -1;
				
				BRec( idx,1+diff*n/2+(n+1)*(i-1) ) = 1;
				BRec( idx,1-diff*n/2+(n+1)*(i+diff-1) ) = -1;
			end
		end
		if mod(diff,2) == 0
			countEq = countEq + n - diff;
		end
	end

	C = (A'*A)\(A'*BTrans);
	%C = inv(A'*A) * A' * BTrans;

	att = C*calPaths(:,1);
	phase = C*calPaths(:,2);
	att = [ -sum(att)
			att ];
	phase = [ -sum(phase)
			  phase ];
	bTransCal = [att phase];
	
	C = (A'*A)\(A'*BRec);

	att = C*calPaths(:,1);
	phase = C*calPaths(:,2);
	att = [ -sum(att)
			att ];
	phase = [ -sum(phase)
			  phase ];
	bRecCal = [att phase];
end
