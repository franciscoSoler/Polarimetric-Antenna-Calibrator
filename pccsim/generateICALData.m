clear all
close all

file2save = 'pruebas/pre';
lenRow = 5000;
fillValue = 128;

% version 1, this scritp generates a QP polarization.
% this only generates one ideal chirp and replicates 

[acq chirpRep]= pccsim_v2( 'row','noise','HV','noOutput' );

passAcqRow = baseb2passb( acq,50e6 );
passChirpRow = baseb2passb( chirpRep,50e6 );

[acq chirpRep]= pccsim_v2( 'panel','noise','HV','noOutput' );

passAcqPanel = baseb2passb( acq,50e6 );
passChirpPanel = baseb2passb( chirpRep,50e6 );


%CE CAL
CE_CIH_RxH = repmat( passChirpRow',1,5 );
CE_CIV_RxV = CE_CIH_RxH;
CE_TxH_CIH_RxH = CE_CIH_RxH;
CE_RxH = CE_CIH_RxH;
CE_TxV_CIV_RxV = CE_CIH_RxH;
CE_RxV = CE_CIH_RxH;

%Noise CAL
RxHV = CE_CIH_RxH;

%ANT CAL
	% ANT CAL Rx
TRM_RxH = [passAcqRow' passAcqPanel'];
CE_CIV_CIH_RxH = CE_CIH_RxH;
TRM_RxV = TRM_RxH;
CE_CIH_CIV_RxV = CE_CIH_RxH;
	
	% ANT CAL Tx
TRM_TxH = TRM_RxH;
CE_TxH_CIH_CIV_RxH = CE_CIH_RxH;
TRM_TxV = TRM_RxH;
CE_TxV_CIV_CIH_RxV = CE_CIH_RxH;


% save data
%PRE cal
%CE CAL
saveSimulatedData( file2save,'w',CE_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIV_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxH_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxV_CIV_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_RxV,lenRow,fillValue );
% Noise CAL
saveSimulatedData( file2save,'a',RxHV,lenRow,fillValue );
%ANT CAL
	% ANT CAL Rx
saveSimulatedData( file2save,'a',TRM_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIV_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',TRM_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIH_CIV_RxV,lenRow,fillValue );
	% ANT CAL Tx
saveSimulatedData( file2save,'a',TRM_TxH,lenRow,fillValue  );
saveSimulatedData( file2save,'a',CE_TxH_CIH_CIV_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',TRM_TxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxV_CIV_CIH_RxV,lenRow,fillValue  );

% post cal
file2save = 'pruebas/post';
%ANT CAL
	% ANT CAL Tx
saveSimulatedData( file2save,'w',TRM_TxH,lenRow,fillValue  );
saveSimulatedData( file2save,'a',CE_TxH_CIH_CIV_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',TRM_TxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxV_CIV_CIH_RxV,lenRow,fillValue  );
	% ANT CAL Rx
saveSimulatedData( file2save,'a',TRM_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIV_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',TRM_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIH_CIV_RxV,lenRow,fillValue );
% Noise CAL
saveSimulatedData( file2save,'a',RxHV,lenRow,fillValue );
%CE CAL
saveSimulatedData( file2save,'a',CE_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_CIV_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxH_CIH_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_RxH,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_TxV_CIV_RxV,lenRow,fillValue );
saveSimulatedData( file2save,'a',CE_RxV,lenRow,fillValue );
