function b = characterizeAntenna( power,qttyAnt )
% this function returns the atenuation of every calibration path, each path
% only uses one pair of transmission / reception TRM.
% input:
%	* power: corresponds to the input power in db and phase, at this stage is unnused

	[x,y] = readCharacterization( qttyAnt );
	b = calculateBi( x,y,qttyAnt );
end
