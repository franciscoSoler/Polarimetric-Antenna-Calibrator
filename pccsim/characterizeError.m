function [ meanErr,stdErr,MSE ] = characterizeError( err,expectedValue )
	meanErr = mean( err );
	stdErr = std( err );
	MSE = sqrt( (expectedValue - meanErr)^2 + stdErr^2 );
end
