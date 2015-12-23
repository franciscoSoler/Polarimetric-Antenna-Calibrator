function [attm, angm] = pccdecoder(acq, chirpRep, walsh_phi_m)

% walsh_phi_m  20x32  : secuencia walsh
% acq          32x2560: 32 lineas de 2560 samples correspondientes al pcc
% chirprep      1x2560: chirp de referencia de loop de la ce 

dec = 1e11;
%%%% DECODIFICACION DE DATOS CRUDOS
rad2deg = 180/pi;
N_elems = 20;   % 7 para paneles, 20 para filas, 140 para mtrs (no se simula promediado por panel/fila)
k = ceil(log2(N_elems));    % Orden de las secuencias de Walsh
N = 2^k;                    % Cantidad de secuencias (y pulsos del modo)
n_lines = N;

%{
signal = repmat( (acq * chirpRep').',N_elems,1 );
integral = signal .* exp( -1i*walsh_phi_m( 1:N_elems,: ));
	% integro todos los tériminos 
signalEst = integral * ones(N,1) /(N*qttyRepeatedRows * (chirpRep*chirpRep'));
	% La integral de arriba tiene dividiendo 2 valores a saber:
	%	N: es ||cj||², TODO: esto está mal tambien, uno tiene que calcular
	%	la correlación de las matrices de generación, la que tiene errores 
	%	con la que no tiene errores.
	%	chirpRep*chirpRep' es el módulo de las dos chirps a la que fue 
	% 	multiplicada la señal, no va la ideal porque no se conoce.

angm = mod( round( angle( signalEst )*dec )/dec*rad2deg,360 ); % signalEst de fase en deg
attm = round( mag2db( abs( signalEst ) )*dec )/dec; % signalEsts de atenuacion en dB power


%line_samples = 4096; %2560
%}
line_samples = length(chirpRep);

for el = 1:N_elems
    m = zeros(n_lines,line_samples);
    for lin = 1:n_lines
        cj = exp(-1i*walsh_phi_m(el,lin)) * conj(chirpRep);  % walsh y replica conjugados
        m(lin,:) = pb2bb(acq(lin,:)) .* cj;  % producto para la correlacion
    end;
    rdo = sum(sum(m))/N;   % integral
    angm(el) = mod(angle(rdo)*rad2deg,360); % rdos de fase en deg
    attm(el) = 20*log10(abs(rdo)/line_samples); % rdos de atenuacion en dB power
end
%}
end