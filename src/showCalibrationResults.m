function showCalibrationResults( idealGainPhaseTrans,realGainPhaseTrans,idealGainPhaseRec,realGainPhaseRec )
   	qttyElements = size( realGainPhaseTrans,1 ) + 2; 
    
    pathsTrans = getAttenuationPath( qttyElements,true );
	estPathTrans = pathsTrans + realGainPhaseTrans;
	idealPathTrans = pathsTrans + idealGainPhaseTrans;
	pathsRec = getAttenuationPath( qttyElements,false );
	estPathRec = pathsRec + realGainPhaseRec;
	idealPathRec = pathsRec + idealGainPhaseRec;
    fprintf('Transmission\n')
    for el = 1:size(pathsTrans, 1)
        fprintf( sprintf( 'Elemento: %2d - Att. est/real: %03.2f/%03.2f - Fase est/real: %06.2f/%06.2f\n',el,estPathTrans( el,1 ),pathsTrans( el,1 ),estPathTrans( el,2 ),pathsTrans( el,2 ) ) );
    end
    fprintf('Reception\n')
    for el = 1:size(pathsRec, 1)
        fprintf( sprintf( 'Elemento: %2d - Att. est/real: %03.2f/%03.2f - Fase est/real: %06.2f/%06.2f\n',el,estPathRec( el,1 ),pathsRec( el,1 ),estPathRec( el,2 ),pathsRec( el,2 ) ) );
    end
    fprintf('Transmission\n')
    for el = 1:size(pathsTrans, 1)
        fprintf( sprintf( 'Elemento: %2d - diffAtt. ideal/real: %03.2f - diffFase ideal/real: %06.2f\n',el,estPathTrans( el,1 ) - idealPathTrans( el,1 ),estPathTrans( el,2 ) - idealPathTrans( el,2 ) ) );
    end
    fprintf('Reception\n')
    for el = 1:size(pathsRec, 1)
        fprintf( sprintf( 'Elemento: %2d - diffAtt. ideal/real: %03.2f - diffFase ideal/real: %06.2f\n',el,estPathRec( el,1 ) - idealPathRec( el,1 ),estPathRec( el,2 ) - idealPathRec( el,2 ) ) );
    end
	x = (1:size(pathsTrans,1))';
	% Grafico de atenuaciones
	figure
	
	subplot( 2,2,1 ); 
	plot(x, estPathTrans( :,1 ),  'Color','r', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, pathsTrans( :,1 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
	legend('Location','Best');
	ylabel('Att [dB]'); xlabel('RMs'); 
	title('Transmission paths');
	grid on

	% Grafico de fases
	subplot( 2,2,2 ); 
	plot(x, estPathTrans( :,2 ), 'Color','r', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, pathsTrans( :,2 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
	legend('Location','Best');
	ylabel('Pha [deg]'); xlabel('RMs'); 
	title('Transmission paths');
	grid on

	subplot( 2,2,3 ); 
	plot(x, estPathRec( :,1 ),  'Color','r', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, pathsRec( :,1 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
	legend('Location','Best');
	ylabel('Att [dB]'); xlabel('RMs'); 
	title('Reception paths');
	grid on

	% Grafico de fases
	subplot( 2,2,4 ); 
	plot(x, estPathRec( :,2 ), 'Color','r', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, pathsRec( :,2 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Real');
	legend('Location','Best');
	ylabel('Pha [deg]'); xlabel('RMs'); 
	title('Reception paths');
	grid on
	
	% la siguiente figura es para las diferencias:
	deltaAtt = 0.5/2;
	deltaPhase = 5.625/2;

	figure
	
	subplot( 2,2,1 ); 
	plot(x, estPathTrans( :,1 ),  'Color','b', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, idealPathTrans( :,1 ), 'Color','g', 'Marker', '.', 'DisplayName', 'Ideal');
	plot(x, idealPathTrans( :,1 )+deltaAtt,  'Color','r', 'Marker', '.', 'DisplayName', 'Upper Limit'); 
	plot(x, idealPathTrans( :,1 )-deltaAtt,  'Color','r', 'Marker', '.', 'DisplayName', 'Lower Limit'); 
	legend('Location','Best');
	ylabel('Att [dB]'); xlabel('RMs'); 
	title('Differences between real/ideal calibrated paths, Transmission');
	grid on

	% Grafico de fases
	subplot( 2,2,2 ); 
	plot(x, estPathTrans( :,2 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, idealPathTrans( :,2 ), 'Color','g', 'Marker', '.', 'DisplayName', 'Ideal');
	plot(x, idealPathTrans( :,2 )+deltaPhase,  'Color','r', 'Marker', '.', 'DisplayName', 'Upper Limit'); 
	plot(x, idealPathTrans( :,2 )-deltaPhase,  'Color','r', 'Marker', '.', 'DisplayName', 'Lower Limit'); 
	legend('Location','Best');
	ylabel('Pha [deg]'); xlabel('RMs'); 
	title('Differences between real/ideal calibrated paths, Transmission');
	grid on

	subplot( 2,2,3 ); 
	plot(x, estPathRec( :,1 ),  'Color','b', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, idealPathRec( :,1 ), 'Color','g', 'Marker', '.', 'DisplayName', 'Ideal');
	plot(x, idealPathRec( :,1 )+deltaAtt,  'Color','r', 'Marker', '.', 'DisplayName', 'Upper Limit'); 
	plot(x, idealPathRec( :,1 )-deltaAtt,  'Color','r', 'Marker', '.', 'DisplayName', 'Lower Limit'); 
	legend('Location','Best');
	ylabel('Att [dB]'); xlabel('RMs'); 
	title('Differences between real/ideal calibrated paths, Reception');
	grid on

	% Grafico de fases
	subplot( 2,2,4 ); 
	plot(x, estPathRec( :,2 ), 'Color','b', 'Marker', '.', 'DisplayName', 'Calibrated'); 
	hold on
	plot(x, idealPathRec( :,2 ), 'Color','g', 'Marker', '.', 'DisplayName', 'Ideal');
	plot(x, idealPathRec( :,2 )+deltaPhase,  'Color','r', 'Marker', '.', 'DisplayName', 'Upper Limit'); 
	plot(x, idealPathRec( :,2 )-deltaPhase,  'Color','r', 'Marker', '.', 'DisplayName', 'Lower Limit'); 
	legend('Location','Best');
	ylabel('Pha [deg]'); xlabel('RMs'); 
	title('Differences between real/ideal calibrated paths, Reception');
	grid on

end

function paths = getAttenuationPath( qttyElements,isTransmission )
    if isTransmission
        [X( :,1 ) X( :,2 )] = textread( 'RFDNcharacterization', '%f %f', qttyElements,'headerlines',0 );
    else
        [X( :,1 ) X( :,2 )] = textread( 'RFDNcharacterization', '%f %f', qttyElements, 'headerlines',qttyElements );
        X = [X( end-1:end,: ); X( 1:end-2,: )];
    end
    paths = [X( 3:end,1 ) + X( 1 ) + X( 2 ) X( 3:end,2 ) + X( 1,2 ) + X( 2,2 )];
end
