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

% Release notes:
% V1:
% - Version original
%
% V2:
% - Misma funcionalidad que V1, con codigo mejor ordenado y comentado


%%%%%%%%%%%%%%%%%%%%%%%%
% INICIALIZACION
%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
deg2rad = pi/180;
rad2deg = 180/pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURACION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cantidad de elementos a caracterizar
N_elems = 140;   % 7 para paneles, 20 para filas, 140 para mtrs (no se simula promediado por panel/fila)

% Fase y atenuacion de cada uno de los N_elems (loop entero sobre la antena=rfdn ida + TRM + rfdn vuelta)
ph_deg = mod((1:N_elems) * 5.625,360); % fase [deg] escalonada, subiendo 5.625deg en cada loop, normalizada de 0 a 360deg
att_db = zeros(1,N_elems);   % atenuacion nula en todos los elementos
att_db(1) = 0.5;           % se setea atenuacion 0.5dB en elemento 1 para poder ver la diferencia en el rdo
att_db(2) = 0.5;
att_db(3) = 0.5;
att_db(4) = 0.5;
att_db(5) = 0.5;
att_db(6) = 0.5;
att_db(7) = 0.5;
att_db(8) = 0.5;
att_db(9) = 0.5;
att_db(10) = 0.5;
att_db(11) = 0.5;
att_db(12) = 0.5;           % se setea atenuacion 1.5dB en elemento 2 para poder ver la diferencia en el rdo
att_db(13) = 0.5;           % se setea atenuacion 2.5dB en elemento 3 para poder ver la diferencia en el rdo
att_db(14) = 0.5;
att_db(15) = 0.5;
att_db(16) = 0.5;
att_db(17) = 0.5;
att_db(18) = 2.5;
att_db(19) = 0.5;
att_db(20) = 0.5;
% Datos del chirp
bw = 10e6;      % ancho de banda [Hz]
fc = 0;    		% frecuencia central [Hz] 0 ->baseband
tp = 20.00e-6;  % duracion [sec]

% Datos de ventana de muestreo
swst = 0e-6;    % retardo de apertura de ventana de muestreo (sampling window start time) [sec]
swl = 2*tp;     % duracion de la ventana de muestreo (sampling window length) [sec]
fs = 120e6;     % frecuencia de sampleo [sec]

%%%% ERRORES EN CHIRP (long loop)

% Setear en 1 para activar error residual de LUT del desfasaje introducido por los TRMs, 0 para desfasaje perfecto
trm_ph_err_on = 0;

% Desviacion estandar del error residual de LUT del desfasaje introducido por los TRMs 
% (importa solo si trm_ph_err_on == 1)
err_ph = 5; % std del error de fase en deg

% Setear en 1 para activar errores en chirp replica
chirp_rep_err_on = 1;

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
icatt = 0.2;  % ical loop path power attenuation [dB] (req. to be known to 0.1dB)
icph = -5;     % ical loop path insertion phase [deg] (req. to be known to 5deg)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA PRINCIPAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% GENERACION DE SECUENCIAS ORTOGONALES PARA LA CODIFICACION PCC
k = ceil(log2(N_elems));    % Orden de las secuencias de Walsh
N = 2^k;                    % Cantidad de secuencias (y pulsos del modo)
% Matriz con las secuencias de desfasajes para cada elemento (filas) 
% en cada pulso (cols) en radianes:
walsh_phi_m = walsh_mtx(k)*pi/2;   
if trm_ph_err_on == 1
 % Se genera una matriz con error de los desfasadores de los TRMs y se agrega a la ideal
 walsh_phi_m_err = walsh_phi_m + randn(size(walsh_phi_m))*err_ph*deg2rad; 
else
 walsh_phi_m_err = walsh_phi_m; % caso ideal sin error
end



% CONSTRUCCION DEL CHIRP QUE RECORRE CADA TRM DE LA ANTENA (ICAL LONG LOOP)
kr = bw/tp;                          % Chirp rate en rango
line_samples = floor(swl*fs);        % Cantidad de samples por linea
t = 0:(1/fs):(line_samples-1)/fs;    % Tiempo para cada sample
f0 = fc-bw/2;                        % Frecuencia inicial
chirp = exp(i*2*pi*(f0*t + kr/2*t.^2)); % chirp ideal, con envolvente unitaria
% Se creo un chirp ideal con envolvente unitaria
% En prox version se puede construir el chirp con errores lineales/cuadraticos/random 
% de la CE(sin hardware de ical)+antena(rfdn ida + trm + rfdn vuelta)



% CONSTRUCCION DEL CHIRP REPLICA DE REFERENCIA (OBTENIDO DE ICAL SHORT LOOP)
icatt = power(10,icatt/-20); % ical loop path attenuation 
icph = icph * deg2rad; % ical loop path phase
if chirp_rep_err_on == 1
 chirprep = icatt * exp(i*icph) * exp(i*2*pi*(f0*t + kr/2*t.^2)); 
else
 chirprep = chirp; % caso con cadena de calibracion CE ideal
end
% Se creo un chirp replica, asumiendo que el short loop (hardware de ical de la CE) agrega 
% un offset fijo de fase y amplitud
% En prox version se puede construir el chirp con errores lineales/cuadraticos/random 
% de la CE(con hardware de ical) 


%%% CONVERSION DE UNIDADES DE ATENUACION Y FASE DE CADA ELEMENTO

% vector con las atenuaciones a medir, correspondientes a c/u de los 
% N_elems elementos o paths, pasados a veces
amp = power(10,att_db/-20);  
% vector con las fases a medir, correspondientes a c/u de los 
% N_elems elementos o paths, pasados a radianes
ph_rad = ph_deg * deg2rad;   



%%%% CONSTRUCCION DE DATOS CRUDOS

n_lines = N;
acq = zeros(n_lines,line_samples);%,'int8');
for lin = 1:n_lines
    linea = zeros(1,line_samples,'double');
    % Generacion de la señal composite correspondiente a cada eco
    for el = 1:N_elems
        % Fase agregada por cada camino (seteo real + codigo walsh agregado, con error del desfasador)
        phi0 = ph_rad(el) + walsh_phi_m_err(el,lin);
        % Construccion de la señal loopeada por el elemento actual
        loop = amp(el) * exp(i*phi0) * chirp;
          % Aca en prox version se podrian agregar señales de error debido a
          % ruido, repetibilidad de switches, isolation en switches de TRMs y power splitters,
          % inestabilidades en hardware de calibracion y nominal
          % (sus reflexiones que generan error en la medicion, porque el camino 
          % nominal y su estabilidad es lo que se mide), y RFI
        % Suma de la señal loopeada por el elemento actual y los otros
        linea = linea +  loop;
    end;
    acq(lin,:)=linea;
end;    

% aca se puede grabar la adquisicion en disco para simularle datos de entrada al 
% procesador de la Processing Chain y validarlo



%%%% DECODIFICACION DE DATOS CRUDOS

for el = 1:N_elems
    m = zeros(n_lines,line_samples);
    for lin = 1:n_lines
        cj = exp(-i*walsh_phi_m(el,lin)) * conj(chirprep);  % walsh y replica conjugados
        m(lin,:) = acq(lin,:) .* cj;  % producto para la correlacion
    end;
    rdo = sum(sum(m))/N;   % integral
    rdo1(el) = rdo/(chirprep*chirprep');
	rdo2(el) = rdo;
	angm(el) = mod(angle(rdo)*rad2deg,360); % rdos de fase en deg
    errp(el) = angm(el) - ph_rad(el)*rad2deg; % error de estimacion de fase
    errp = (errp>180)*(-360)+errp; % se lleva todo alrededor de 0
    errp = (errp<-180)*(+360)+errp; % se lleva todo alrededor de 0
    attm(el) = -20*log10(abs(rdo)/line_samples); % rdos de atenuacion en dB power    
    erra(el) = attm(el) - att_db(el);    
    %erra(el) = -20*log10( ampm(el) -  amp(el) );
end;   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT DE RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Atenuacion/fase estimada vs real
for el = 1:N_elems
    disp(sprintf('Elemento: %2d - Att. est/real: %03.2f/%03.2f - Fase est/real: %06.2f/%06.2f',el,attm(el),att_db(el),angm(el),ph_deg(el)))
end;  

% Error de estimacion de atenuacion/fase
for el = 1:N_elems
    disp(sprintf('Elemento: %2d - Error fase/amp: %05.4fdeg/%05.4fdB',el,errp(el),erra(el)));    
end;

% Estadisticas de error
mp = mean(errp); 
sp = std(errp); 
rssp = sqrt(mp^2+sp^2);
ma = mean(erra); 
sa = std(erra); 
rssa = sqrt(ma^2+sa^2);
disp(sprintf('Pha estimation error mean/std/rss: %02.2f / %02.2f / %02.2f deg',mp,sp,rssp))
disp(sprintf('Att estimation error mean/std/rss: %02.2f / %02.2f / %02.2f dB ',ma,sa,rssa))

% Inicializacion de vectores auxiliares para graficos
x = 1:N_elems;
y = ones(1,N_elems);
close all

% Grafico de atenuaciones
figure; 
subplot(4,1,[1 2]); hold on;  
plot(x, attm,  'Color','r', 'Marker', '.', 'DisplayName', 'Estimated'); 
plot(x, att_db, 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
legend('Location','Best');
ylabel('Att [dB]'); xlabel('Panel,Row or TRM level path'); 
subplot(4,1,[3 4]); hold on; 
plot(x, erra, 'Color','r', 'Marker', '.', 'DisplayName', 'Difference/path');  
plot(x, y*ma, 'Color','g', 'DisplayName', sprintf('Dif. media = %02.2fdB',ma));  
plot(x, y*sa, 'Color','b', 'DisplayName', sprintf('Dif. std = %02.2fdB',sa));  
legend('Location','Best');
xlabel('Panel,Row or TRM level path'); 
ylabel(sprintf('Att [dB] \n Error sqrt(media^2+std^2)=%02.2fdB',rssa))

% Grafico de fases
figure; 
subplot(4,1,[1 2]); hold on; 
plot(x, angm,  'Color','r', 'Marker', '.', 'DisplayName', 'Estimated'); 
plot(x, ph_deg, 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
legend('Location','Best');
ylabel('Pha [deg]'); xlabel('Panel,Row or TRM level path'); 
subplot(4,1,[3 4]); hold on; 
plot(x, errp, 'Color','r', 'Marker', '.', 'DisplayName', 'Difference/path');  
plot(x, y*mp, 'Color','g', 'DisplayName', sprintf('Dif. media = %02.2fdeg',mp));  
plot(x, y*sp, 'Color','b', 'DisplayName', sprintf('Dif. std = %02.2fdeg',sp));  
legend('Location','Best');
xlabel('Panel,Row or TRM level path'); 
ylabel(sprintf('Pha [deg] \n Error sqrt(media^2+std^2)=%02.2fdeg',rssp))


