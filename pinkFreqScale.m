function [scaledFreqs] = pinkFreqScale(freqMin,freqMax,numFreqs,exponent)

% 1/f equation:
% pwr = 1 ./ (f .^ exponent)

% Recommended exponent is 0.05 based on Voytek & Gazzaley 2015

% find function value at min and max frequency of interest
% (max power at min freq and vice versa)
maxPower = 1 ./ (freqMin .^ exponent);
minPower = 1 ./ (freqMax .^ exponent);

% linear spacing in log-space
powerSteps = linspace(maxPower,minPower,numFreqs);

% solve the inverse of the 1/f equation to find the scaled frequency
scaledFreqs = nthroot((1 ./ powerSteps),exponent);

end % end function