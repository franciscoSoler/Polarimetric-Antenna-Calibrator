function [walkn] = walsh_seq(k,n)
%[walkn] = walsh_seq(k,n)
%  walkn    Enesima secuencia de Walsh de orden k, dimension (1,2^k)
%  k        Orden de la secuencia de Walsh. k>0
%  n        Enesima secuencia. n>0

	N = 2^k;

	if n == 0
    	walkn = ones(1,N);    
		return
	end
	m = 1 + floor(log2(n)); 
	ii = 0:N-1;
	rad = (-1).^floor( ii* 2^m /N );
	walkn = rad .* walsh_seq(k,2^m-1-n);
end
