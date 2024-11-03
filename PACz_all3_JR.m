function [PACz_all3,PAC_all3,iteratedPAC_all3,ampByPhase] = PACz_all3_JR(lowFreq_phase,highFreq_power,highFreq_lowPhase,varargin)
%% Run phase amplitude coupling analysis
% Written by Justin Riddle, PhD
%   Inquiries to justin_riddle@med.unc.edu
%
% Calculates phase amplitude coupling according to three different methods
% and calculates normalized value from null-distribution created by
% randomly shifting the input data
% INPUTS:
%   lowFreq_phase = time series of phase values of low frequency
%   highFreq_power = time series of amplitude values of high frequency
%   highFreq_lowPhase = phase values of low frequency conv with high amp
% 
% OPTIONAL INPUTS:
%   options = structure with possible fields:
%       numPhaseBins = double, number of bins for average amp in phase
%       numIterations = double, change iteraction number for null dist
%
% Methods description:
%   1. Mike X. Cohen method from Analyzing Neural Time Series Data
%   2. Tort calculation
%   3. Ryan Canolty method

METHODS = {'cohen','tort','canolty'};
numMethods = length(METHODS);

%% options that can be tweaked
if ~isempty(varargin)
    options = varargin{1};
else
    options = struct([]);
end
if isfield(options,'numPhaseBins')
    numPhaseBins = options.numPhaseBins;
else
    numPhaseBins = 30;
end
if isfield(options,'numIterations')
    numIterations = options.numIterations;
else
    numIterations = 1000;
end

%% COHEN VERSION OF PAC
% calculate PAC value for this subject,region,pacpair% M*e^(i*theta) polar notation for a complex signal
    % The PAC formula here takes the phase of the low
    % frequency wavelet and amplitude of the high
    % frequency wavelet - creates a new signal that is
    % those two combined. Then take the mean across a
    % window of time of this new hybrid signal
    % after mean you have a complex number, now
    % calculate the magnitude of that
PAC_cohen = abs(mean(highFreq_power.*exp(1i*lowFreq_phase)));

%% power by phase
ampByPhase = NaN(1,numPhaseBins);
phaseSpacing=linspace(min(lowFreq_phase),max(lowFreq_phase),numPhaseBins+1);
for binIdx=1:numPhaseBins
    ampByPhase(binIdx) = mean(highFreq_power(lowFreq_phase>phaseSpacing(binIdx) & lowFreq_phase<phaseSpacing(binIdx+1)));
end

%% TORT VERSION OF PAC
% Calculate entropy of average amplitude in the phase bins
pAmp = ampByPhase./sum(ampByPhase);
entropyAmp = -1 * sum(pAmp .* log(pAmp));
PAC_tort = 1 - (entropyAmp./log(numPhaseBins));

%% CANOLTY VERSION OF PAC
angleDiff = lowFreq_phase - highFreq_lowPhase;
% Phase locking value
PAC_canolty = abs(nanmean(exp(1i*angleDiff)));

% All three version are calculated
PAC_all3  = [PAC_cohen PAC_tort PAC_canolty];

% Running iterative boot-strapping
iteratedPAC_all3 = NaN(numIterations,numMethods);
numTotalTimePoints = length(lowFreq_phase);

% must move the data at least 10% (max 90%) for boot-strapping
shift = 0.1;
% shift the power values from 10% to 90% for each
% iteration, then calculate PAC
minShiftNum = ceil(shift*numTotalTimePoints);
numPossRand = floor(numTotalTimePoints*(1-(shift+shift)));

for iterIdx = 1:numIterations
    
    %% Randomly shift the data
    % randomly select a number of time points and move the data
    randIdx = randi(numPossRand);
    shiftPoints = minShiftNum + randIdx;
    
    % Shift high frequency power
    thisIter_highFreq_pwr = [highFreq_power(shiftPoints:end) ...
        highFreq_power(1:(shiftPoints-1))];

    % Shift high frequency low phase
    thisIter_highFreq_lowPhase = [highFreq_lowPhase(shiftPoints:end) ...
        highFreq_lowPhase(1:(shiftPoints-1))];
    
    
    %% Calculate Cohen's PAC
    iter_PAC_cohen = abs(mean(thisIter_highFreq_pwr.*exp(1i*lowFreq_phase)));
    iteratedPAC_all3(iterIdx,1) = iter_PAC_cohen;
    
    
    %% Calculate Tort's PAC
    iter_ampByPhase = NaN(1,numPhaseBins);
    for binIdx=1:numPhaseBins
        iter_ampByPhase(binIdx) = mean(thisIter_highFreq_pwr(...
            lowFreq_phase>phaseSpacing(binIdx)...
            & lowFreq_phase<phaseSpacing(binIdx+1)));
    end
    pAmp = iter_ampByPhase./sum(iter_ampByPhase);
    entropyAmp = -1 * sum(pAmp .* log(pAmp));
    iter_PAC_tort = 1 - (entropyAmp./log(numPhaseBins));
    iteratedPAC_all3(iterIdx,2) = iter_PAC_tort;
    
    
    %% Calcualte Canolty's PAC (sometimes called Modulation Index)
    angleDiff = lowFreq_phase - thisIter_highFreq_lowPhase;
    % Phase locking value
    iter_PAC_canolty = abs(nanmean(exp(1i*angleDiff)));
    iteratedPAC_all3(iterIdx,3) = iter_PAC_canolty;
end

%% Normalize the data
% Z-score the PAC value based on the iterated values
PACz_cohen   = (PAC_cohen   - mean(iteratedPAC_all3(:,1))) / std(iteratedPAC_all3(:,1));
PACz_tort    = (log(PAC_tort)    - mean(log(iteratedPAC_all3(:,2)))) / std(log(iteratedPAC_all3(:,2)));
PACz_canolty = (PAC_canolty - mean(iteratedPAC_all3(:,3))) / std(iteratedPAC_all3(:,3));

%% Store for output
PACz_all3 = [PACz_cohen PACz_tort PACz_canolty];

end