% first, the attenuation of each path is getted
calPaths = characterizeAntenna(0,5);

%[diffTrans diffRec] = test( calPaths );
[diffTrans1 diffRec1] = getSignalShift( calPaths );
[idealGainPhaseTrans, realGainPhaseTrans] = calibrateAntenna( diffTrans );
[idealGainPhaseRec, realGainPhaseRec] = calibrateAntenna( diffRec );
showCalibrationResults( idealGainPhaseTrans,realGainPhaseTrans,idealGainPhaseRec,realGainPhaseRec )
%{
A = [2  0  0  0  1  1  0  0  0 -1  0  0  0
     2 -2  0  0  0  0  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0  0  0  0  0
     0  2  0  0  1  0  1  0  0  0 -1  0  0
     2  0 -2  0  0  0  0  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0  0  0  0  0
     0  0  2  0  1  0  0  1  0  0  0 -1  0
     2  0  0 -2  0  0  0  0  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0  0  0  0  0
     0  0  0  2  1  0  0  0  1  0  0  0 -1
     0  1  0  0  0  1  0  0  0  0  0  0  0
     1  0 -1  0  0  0  0  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0  0  0  0  0
     0  0  0  1  0  0  1  0  0  0  0  0  0
     0  0  2  0  0  1  1  0  0  0  0  0  0
     0  0  0  2  0  1  0  1  0  0  0  0  0];
 % el rango es de 10, necesito 3 dimensiones mas grrgr
rank(A)


%a b c d b1 b2 b3 b4 e f g h b21
B = [1  0  0  0  1  0  0  0  1  0  0  0  0
     2 -2  0  0  0  0  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0  0  0  0  0
     0  1  0  0  1  0  0  0  0  1  0  0  0
     2  0 -2  0  0  0  0  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0  0  0  0  0
     0  0  1  0  1  0  0  0  0  0  1  0  0
     2  0  0 -2  0  0  0  0  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0  0  0  0  0
     0  0  0  2  1  0  0  0  0  0  0  1  0
     1  0 -1  0  0  0  0  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0  0  0  0  0
     0 -1  0  2  0  0  0  1  0  0  0  0  0];
rank(B)
% el rango es de 8
%

%a b c d b1 b4 e f g h
B = [1  0  0  0  1  0  1  0  0  0
     2 -2  0  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0  0
     0  1  0  0  1  0  0  1  0  0
     2  0 -2  0  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0  0
     0  0  1  0  1  0  0  0  1  0
     2  0  0 -2  0  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0  0
     0  0  0  2  1  0  0  0  0  1
     1  0 -1  0  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0  0
     0 -1  0  2  0  1  0  0  0  0];
rank(B)
% el rango es de 8, 2 menos del que necesito

%
Complemento calculos hechos con las 8 hojas
b4 se la calcula como
b4 = G + yd - ya + b1 = G + yd + x1 + x5
como x1 y x5 son iguales y conocidos, la parte de incógnitas son
b4 = G + yd



asumiendo yd e ya conocidos, conozco b1

%a b c d b1 e f g h yd
B = [1  0  0  0  1  1  0  0  0  0
     2 -2  0  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0  0
     0  1  0  0  1  0  1  0  0  0
     2  0 -2  0  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0  0
     0  0  1  0  1  0  0  1  0  0
     2  0  0 -2  0  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0  0
     0  0  0  2  1  0  0  0  1  0
     1  0 -1  0  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0  0
     0 -1  0  2  0  0  0  1  0  1];
rank(B)

% asumiendo deltade y b1 conocido:

%

%a b c d e f g h
B = [1  0  0  0  1  0  0  0
     2 -2  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0
     0  0  2 -2  0  0  0  0
     0  1  0  0  0  1  0  0
     2  0 -2  0  0  0  0  0
     0  2  0 -2  0  0  0  0
     0  0  1  0  0  0  1  0
     2  0  0 -2  0  0  0  0
     2  0  0 -2  0  0  0  0
     0  0  0  2  0  0  0  1
     1  0 -1  0  0  0  0  0
     0  1  0 -1  0  0  0  0
     0 -1  2 -1  0  0  0  0
     0 -1  0  2  0  0  1  0];
rank(B)

% NOT SUCCESS
%

% ahora digo que d y h eran conocidos
%a b c e f g 
B = [1  0  0  1  0  0
     2 -2  0  0  0  0
     0  2 -2  0  0  0
     0  0  2  0  0  0
     0  1  0  0  1  0
     2  0 -2  0  0  0
     0  2  0  0  0  0
     0  0  1  0  0  1
     2  0  0  0  0  0
     2  0  0  0  0  0
     0  0  0  0  0  0
     1  0 -1  0  0  0
     0  1  0  0  0  0
     0 -1  2  0  0  0
     0 -1  0  0  0  1];
rank(B)
% seguro este último paso está mal, tengo que justificarlo con ecuaciones
% pero estoy cansado y esta blaaa

%
% partiendo de la tabla B de la línea 79, si digo que 

b1 = e + f + g + h - yd
%

%a b c d b1 e f g h yd
B = [1  0  0  0  1  1  0  0  0  0
     2 -2  0  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0  0
     0  1  0  0  1  0  1  0  0  0
     2  0 -2  0  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0  0
     0  0  1  0  1  0  0  1  0  0
     2  0  0 -2  0  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0  0
     0  0  0  2  1  0  0  0  1  0
     1  0 -1  0  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0  0
     0 -1  0  2  0  0  0  1  0  1
	 0  0  0  0 -1  1  1  1  1 -1];
rank(B)

% paso de rango 8 aaaaa rango 9 :D
% si en vez de agregar la ecuación simplemente la reemplazo


%a b c d e f g h yd
B = [1  0  0  0  2  1  1  1 -1
     2 -2  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0
     0  1  0  0  1  2  1  1 -1
     2  0 -2  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0
     0  0  1  0  1  1  2  1 -1
     2  0  0 -2  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0
     0  0  0  1  1  1  1  2 -1
     1  0 -1  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0
     0 -1  0  2  0  0  1  0  1];
rank(B)
% rango = 8, falta 1 ecuación
% llego a lo mismo pero eliminando la ecuacion, cual es la siguiente a elimiar??

%a b c d e f g h yd
B = [1  0  0  0  2  1  1  1
     2 -2  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0
     0  0  2 -2  0  0  0  0
     0  1  0  0  1  2  1  1
     2  0 -2  0  0  0  0  0
     0  2  0 -2  0  0  0  0
     0  0  1  0  1  1  2  1
     2  0  0 -2  0  0  0  0
     2  0  0 -2  0  0  0  0
     0  0  0  1  1  1  1  2
     1  0 -1  0  0  0  0  0
     0  1  0 -1  0  0  0  0
     0 -1  2 -1  0  0  0  0
     0 -1  0  2  0  0  1  0];
rank(B)
% si digo que yd es conocid llego a la solucion :D
% aca repito pero eliminando yd con la ecuacion extra
%a b c d b1 e f g h
B = [1  0  0  0  1  1  0  0  0
     2 -2  0  0  0  0  0  0  0
     0  2 -2  0  0  0  0  0  0
     0  0  2 -2  0  0  0  0  0
     0  1  0  0  1  0  1  0  0
     2  0 -2  0  0  0  0  0  0
     0  2  0 -2  0  0  0  0  0
     0  0  1  0  1  0  0  1  0
     2  0  0 -2  0  0  0  0  0
     2  0  0 -2  0  0  0  0  0
     0  0  0  1  1  0  0  0  1
     1  0 -1  0  0  0  0  0  0
     0  1  0 -1  0  0  0  0  0
     0 -1  2 -1  0  0  0  0  0
     0 -1  0  2  0  0  0  1  0
	 0  0  0  0 -1  1  1  1  1];
rank(B)

% SUCCESS, NO SE QUE ONDA CON LO DE B1, tengo que probar si funciona este 
% método


% Trato de buscar la que es LD
B = [1  2  1  1  1 -1
	 1  1  2  1  1 -1
	 1  1  1  2  1 -1
	 1  1  1  1  2 -1
	 1  0  0  1  0  1];

rank(B)

%}

