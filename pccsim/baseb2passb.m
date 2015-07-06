function passBand = baseb2passb( signalbb,fc )
% this function convert the baseband signal to passband signal
% inputs:
%	signalbb: is the baseband signal
%	fc: is the passBand central frequency
	XA = signalbb*exp(1i*2*pi*fc);
	passBand = real(XA);
end
