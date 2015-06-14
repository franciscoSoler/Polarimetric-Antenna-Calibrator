function [idealTRMValue, realTRMValue] = calibrateAntenna( pathDiff )
% pathDiff corresponds to attenuation and phase shift of each transmission 
% or reception path into the antenna
% this function is dummie, are noted that the attenuators and phase shifter
% have their limits. Could be the case of the method tries to go lower than 
% the attenuator could go, the same of the phase shifter.
% limits:
%	attenuator = 0 - 31.5 dB (Amplitude settings range) 
%	phase shift = it haven't (Phase settings range 360º)
% steps:
% 	attenuator = 0.5 dB
%	phase shift = 5.625º
	attStep = 0.5;
	phaseStep = 5.625;
    realTRMValue = [-round( pathDiff( :,1 )/attStep )*attStep -round( pathDiff( :,2)/phaseStep )*phaseStep];
    idealTRMValue = - pathDiff

%	en transmision el atenuador se FIJA para tener los 4 dBcp
end
