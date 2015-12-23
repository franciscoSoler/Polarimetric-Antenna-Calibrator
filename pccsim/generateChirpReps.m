%this script only generates the chirpRep for each polarization
close all
clear all

file2save = 'pruebas/chirps';
lenRow = 5000;
fillValue = 128;


[acq chirpRep]= pccsim_v2( 'row','noise','HV','noOutput' );

passChirpRow = baseb2passb( chirpRep,50e6 );
chirp2save = repmat( passChirpRow',1,4 );


saveSimulatedData( file2save,'w',chirp2save,lenRow,fillValue );
