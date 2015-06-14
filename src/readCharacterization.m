function [X, Y] = readCharacterization( qttyAnt )
	% in this system in particular
	qttyLines = 2*qttyAnt + 4;
	[X(:,1) X(:,2)] = textread('RFDNcharacterization', '%f %f', qttyLines);
	[Y(:,1) Y(:,2)] = textread('RMcharacterization', '%f %f', qttyAnt);
end
