function [acq chirpRep phi0 amp] = pccsim_v2( level,mode,pol,showOutput ) 
%%%%%%%%%%%%%%%%%%%%%%%%
% Simulador de PCC
%%%%%%%%%%%%%%%%%%%%%%%%

% Descripcion:
% Implementacion de simulador de PCC en base a la publicacion de DLR y
% Austrium "Individual TR Module Characterisation of the TerraSAR-X 
% antenna by calibration pulse sequences with orthogonal codes"
% El software simula los pulsos generados y codificados que forman la 
% senal composite de PCC, simulando los datos crudos de una adquisicion,
% introduciendo errores en el medio que se pueden activar y desactivar, 
% y despues hace la decodificacion PCC comparando los valores originales
% con los estimados para calcular el error. 
% Inputs:
%	level: it can be row, panel or TRM.
%	mode: Tx, Rx or noise
%	pol: H, V or HV

% Release notes:
% V1:
% - Version original
%
% V2:
% - Misma funcionalidad que V1, con codigo mejor ordenado, comentado y 
%	optimizado


%%%%%%%%%%%%%%%%%%%%%%%%
% INICIALIZACION
%%%%%%%%%%%%%%%%%%%%%%%%
dec = 1e11;
deg2rad = pi/180;
rad2deg = 180/pi;
configFile = 'PCCconfig';

validLvl = strcmp( level,'row' ) || strcmp( level,'panel' ) || ...
			strcmp( level,'TRMs' );
if ~validLvl
	error( 'the level parameter must be one of the next tree: row, panelor TRMs' );
end

validMode = strcmp( mode,'Tx' ) || strcmp( mode,'Rx' ) || ...
			strcmp( mode,'noise' );
if ~validMode
	error( 'the mode parameter must be one of the next tree: Tx, Rx or noise' );
end

validPol = strcmp( pol,'H' ) || strcmp( pol,'V' ) || strcmp( pol,'HV' );
if ~validPol
	error( 'the level parameter must be one of the next tree: H, V or HV' );
end

notValidPolHV = strcmp( pol,'HV') && ~strcmp( mode,'noise' );
if notValidPolHV
	error( 'The pol HV is only used for noise calibration' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURACION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cantidad de elementos a caracterizar (casi casi, falta terminar, no esta del todo bien)
config = loadjson( configFile );

% cantidad de paneles * cantidad de filas (140 TRMs)(no se simula promediado 
% por panel/fila)

N_elems = config.antennaConfig.rows * config.antennaConfig.panels;


% Fase y atenuacion de cada uno de los N_elems (loop entero sobre la antena=rfdn ida + TRM + rfdn vuelta)
ph_deg = mod( ( 1:N_elems )' * config.steeringAngle,360 ); 
% fase [deg]  escalonada, subiendo 5.625deg en cada loop, normalizada de 0 a 
% 360deg
att_db = full( config.attElementsDb );

% Datos del chirp
bw = config.chirpData.bw;  % ancho de banda [Hz]
fc = config.chirpData.fc;  % frecuencia central [Hz] 0 ->baseband
tp = config.chirpData.tp;  % duracion [sec]

% Datos de ventana de muestreo
swst = config.samplingWindow.swst; % retardo de apertura de ventana de 
								% muestreo (sampling window start time) [sec]
fs = config.samplingWindow.fs;		% frecuencia de sampleo [sec]
swl = 2*tp; % duracion de la ventana de muestreo (sampling window length) [sec]

%%%% ERRORES EN CHIRP (long loop)

% Setear en 1 para activar error residual de LUT del desfasaje introducido por los TRMs, 0 para desfasaje perfecto
trm_ph_err_on = config.errors.trmPhErrOn;

% Desviacion estandar del error residual de LUT del desfasaje introducido por los TRMs 
% (importa solo si trm_ph_err_on == 1)
err_ph = config.errors.phError * deg2rad; % std del error de fase en rad

% Setear en 1 para activar errores en chirp replica
chirp_rep_err_on = config.errors.chirpRepErrOn;

% Se asume en esta version solo el error de codificacion, con chirp sin 
% errores lineales/cuadraticos/random
% En prox version se pueden agregar señales de error debido a
% ruido, repetibilidad de switches, isolation en switches de TRMs y power splitters,
% inestabilidades en hardware de calibracion y nominal
% (sus reflexiones que generan error en la medicion, porque el camino 
% nominal y su estabilidad es lo que se mide), y RFI



%%%% ERRORES EN REFERENCE CHIRP REPLICA (short loop)

% Errores del chirp replica obtenido de ICAL short loop con respecto al long loop
% (importa solo si trm_ph_err_on == 1)
% Se asume por esta version un error en amplitud y fase como offsets fijos 
% sin variaciones lineales/cuadrat/random. 
% En esto impactan igual que en lo anterior, ruido, repetibilidad de swithces, 
% isolation, inestabilidades del hardware, error residual de caracterizacion 
% en caso de que este caracterizado, RFI.  
icatt = config.errorsICAL.att;  % ical loop path power attenuation [dB] (req.
								% to be known to 0.1dB)
icph = config.errorsICAL.phShift; % ical loop path insertion phase [deg] (req.
								  % to be known to 5deg)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA PRINCIPAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% GENERACION DE SECUENCIAS ORTOGONALES PARA LA CODIFICACION PCC

if strcmp( level,'row' ) % 20 elementos
	nEl = config.antennaConfig.rows;
	qttyRepeatedRows = config.antennaConfig.panels;
	ph_deg = mod( ( ceil( (1:N_elems)/qttyRepeatedRows ) )' * config.steeringAngle,360 ); 
elseif strcmp( level,'panel' ) 
	nEl = config.antennaConfig.panels;
	qttyRepeatedRows = config.antennaConfig.rows;
	ph_deg = mod( ( ceil( (1:N_elems)/qttyRepeatedRows ) )' * config.steeringAngle,360 ); 
else
	nEl = N_elems;
	qttyRepeatedRows = 1;
	ph_deg = mod( ( 1:N_elems )' * config.steeringAngle,360 ); 
end

k = ceil(log2(nEl));    	% Orden de las secuencias de Walsh
N = 2^k;                    % Cantidad de secuencias (y pulsos del modo)
% Matriz con las secuencias de desfasajes para cada elemento (filas) 
% en cada pulso (cols) en radianes:
walsh_phi_m = createPhWalshMtx( 'order',k,'phShift',pi/2,'duplicateRows',qttyRepeatedRows );   
if trm_ph_err_on
	% Se genera una matriz con error de los desfasadores de los TRMs y se agrega a la ideal
	walsh_phi_m_err = createPhWalshMtx( 'walshPhMtx',walsh_phi_m,'phErr',err_ph );
else
	walsh_phi_m_err = walsh_phi_m; % caso ideal sin error
end



% CONSTRUCCION DEL CHIRP QUE RECORRE CADA TRM DE LA ANTENA (ICAL LONG LOOP)

chirp = createChirp( fs,fc,bw,tp,swl ); % chirp ideal, con envolvente unitaria
% Se creo un chirp ideal con envolvente unitaria
% En prox version se puede construir el chirp con errores lineales/cuadraticos/random 
% de la CE(sin hardware de ical)+antena(rfdn ida + trm + rfdn vuelta)



% CONSTRUCCION DEL CHIRP REPLICA DE REFERENCIA (OBTENIDO DE ICAL SHORT LOOP)
if chirp_rep_err_on
	chirpRep = createChirp( fs,fc,bw,tp,swl,'attenuation',icatt,'phaseShift',icph );
else
	chirpRep = chirp; % caso con cadena de calibracion CE ideal
end
% Se creo un chirp replica, asumiendo que el short loop (hardware de ical de la CE) agrega 
% un offset fijo de fase y amplitud
% En prox version se puede construir el chirp con errores lineales/cuadraticos/random 
% de la CE(con hardware de ical) 


%%% CONVERSION DE UNIDADES DE ATENUACION Y FASE DE CADA ELEMENTO

% vector con las atenuaciones a medir, correspondientes a c/u de los 
% N_elems elementos o paths, pasados a veces
amp = db2mag(-att_db);
% vector con las fases a medir, correspondientes a c/u de los 
% N_elems elementos o paths, pasados a radianes
ph_rad = ph_deg * deg2rad;   



%%%% CONSTRUCCION DE DATOS CRUDOS

% Fase agregada por cada camino (seteo real + codigo walsh agregado, con 
% error del desfasador)
phi0 = repmat( ph_rad,1,N ) + walsh_phi_m_err( 1:N_elems,: );
% Construccion de la señal loopeada por cada elemento y sumada entre todos 
acq = (amp' * exp( 1i*phi0 )).' * chirp;
		% TODO: debería hacer chirp * lo otro, acq me queda traspuesta.
% Aca en prox version se podrian agregar señales de error debido a ruido, 
% repetibilidad de switches, isolation en switches de TRMs y power splitters,
% inestabilidades en hardware de calibracion y nominal
% (sus reflexiones que generan error en la medicion, porque el camino 
% nominal y su estabilidad es lo que se mide), y RFI

% aca se puede grabar la adquisicion en disco para simularle datos de entrada al 
% procesador de la Processing Chain y validarlo

if strcmp( showOutput,'noOutput')
	return
end


%%%% DECODIFICACION DE DATOS CRUDOS

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
	% error de estimacion de fase llevado a 0	
errp = mod( angm - ph_deg + 180,360 ) - 180;
attm = round(-mag2db( abs( signalEst ) )*dec )/dec; % signalEsts de atenuacion en dB power

erra = attm - att_db;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT DE RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Atenuacion/fase estimada vs real
for el = 1:N_elems
    disp( sprintf( 'Elemento: %2d - Att. est/real: %03.2f/%03.2f - Fase est/real: %06.2f/%06.2f',el,attm( el ),att_db( el ),angm( el ),ph_deg( el ) ) )
end  

% Error de estimacion de atenuacion/fase
for el = 1:N_elems
    disp( sprintf( 'Elemento: %2d - Error fase/amp: %05.4fdeg/%05.4fdB',el,errp( el ),erra( el ) ) );    
end

% Estadisticas de error
[mp sp rssp] = characterizeError( errp,0 ); 
[ma sa rssa] = characterizeError( erra,0 ); 

disp(sprintf('Pha estimation error mean/std/rss: %02.2f / %02.2f / %02.2f deg',mp,sp,rssp))
disp(sprintf('Att estimation error mean/std/rss: %02.2f / %02.2f / %02.2f dB ',ma,sa,rssa))

% Inicializacion de vectores auxiliares para graficos
x = (1:N_elems)';
y = ones( N_elems,1 );
close all

% Grafico de atenuaciones
figure 
subplot(2,1,1); 
plot(x, attm,  'Color','r', 'Marker', '.', 'DisplayName', 'Estimated'); 
hold on
plot(x, att_db, 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
legend('Location','Best');
ylabel('Att [dB]'); xlabel('Panel,Row or TRM level path'); 

subplot(2,1,2); 
plot(x, erra, 'Color','r', 'Marker', '.', 'DisplayName', 'Difference/path');  
hold on
plot(x, y*ma, 'Color','g', 'DisplayName', sprintf('Dif. media = %02.2fdB',ma));
plot(x, y*sa, 'Color','b', 'DisplayName', sprintf('Dif. std = %02.2fdB',sa));  
legend('Location','Best');
xlabel('Panel,Row or TRM level path'); 
ylabel(sprintf('Att [dB] \n Error sqrt(media^2+std^2)=%02.2fdB',rssa))

% Grafico de fases
figure 
subplot( 2,1,1 ); 
plot(x, angm, 'Color','r', 'Marker', '.', 'DisplayName', 'Estimated'); 
hold on
plot(x, ph_deg, 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
legend('Location','Best');
ylabel('Pha [deg]'); xlabel('Panel,Row or TRM level path'); 

subplot( 2,1,2 ); 
plot(x, errp, 'Color','r', 'Marker', '.', 'DisplayName', 'Difference/path');  
hold on
plot(x, y*mp, 'Color','g', 'DisplayName', sprintf('Dif. media = %02.2fdeg',mp));
plot(x, y*sp, 'Color','b', 'DisplayName', sprintf('Dif. std = %02.2fdeg',sp));
legend('Location','Best');
xlabel('Panel,Row or TRM level path'); 
ylabel(sprintf('Pha [deg] \n Error sqrt(media^2+std^2)=%02.2fdeg',rssp))
