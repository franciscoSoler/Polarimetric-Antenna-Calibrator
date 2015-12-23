function [walkm] = walsh_mtx(k)
%  walkm    MatrizEnesima secuencia de Walsh de orden k, dimension (2^k,2^k)
%  k        Orden de la secuencia de Walsh. k>=0
	walkm = 1;

	for i = 1:k
		walkm = [ walkm  walkm
				  walkm -walkm ];
	end
end
