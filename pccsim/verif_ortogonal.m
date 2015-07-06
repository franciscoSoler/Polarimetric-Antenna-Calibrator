% Verificacion de ortogonalidad de las secuencias de Walsh generadas


k=3;

walsh_m = walsh_mtx(k);

% Verificacion de ortogonalidad
% Se correlaciona cada secuencia (fila de la matriz) con si misma y las
% demas. corre(n,m) es el rdo. de la correlacion (en tau=0) entre la fila n
% y la fila m. 

corre = walsh_m * walsh_m'

walsh_m
walsh_ph_m = exp(i*walsh_m*pi/2)
corre2 = walsh_ph_m * walsh_ph_m'
