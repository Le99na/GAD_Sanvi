function [freqPhase,freqAmp]=phaseAmp_calculator(sampleRate,freq,inputData,numCycles,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CREDITS                                                               %
%                                                                         %
%       By Justin Riddle, PhD  & Sankaraleengam Alagapan, PhD             %
%       Flavio Frohlich's Lab                                             %
%       University of North Carolina at Chapel Hill                       %
%       Last updated on October 20th, 2018                                %
%       Email for contact: justin_riddle@med.unc.edu                      %
%                                                                         %
%% Runs amplitude and phase estimate on time series data                  %
% INPUTS:                                                                 %
%   sampleRate: double for the rate at which the data was sampled         %
%   freq: what is the frequency for this estimate                         %
%   inputData: time series of data to estimate spectral amp and phase     %
%   numCycles: the number of cycles for wavelet convolution               %
%   varargin: optional input of a precomputed wavelet                     %
%       this will save time with many function calls                      %
%                                                                         %
% OUTPUT:                                                                 %
%   freqPhase: phase of the frequency for the time series                 %
%   freqAmp: amplitude of the frequency of the time series                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the user has provided a wavelet to use
if length(varargin)==1
    wavelet = varargin{1};
else % otherwise geneate that wvelet
    numOscPerSide = numCycles /2;
    waveletF = freq;
    waveletT = 1./waveletF;
    waveletS = (2*numOscPerSide)./(2*pi*waveletF);
    timePerSide = numOscPerSide * waveletT;
    waveletTime = -timePerSide: 1/sampleRate : timePerSide;
    waveletCosine = cos( 2*pi*waveletF*waveletTime );
    waveletSine = sin( 2*pi*waveletF*waveletTime );
    waveletGaussian = (waveletS^(-0.5)*pi^(-0.25))* exp(-waveletTime.^2/(2*waveletS^2));
    waveletReal = waveletCosine.*waveletGaussian;
    waveletImag = waveletSine.*waveletGaussian;
    wavelet = waveletReal+1i*waveletImag;
end
% assert an odd nubmer of cycles for the wavelet
assert(mod(numCycles,2)==1);

% Code is partially adapted from Mike X. Cohen's textbook:
%       "Analyzing Neural Time Series Data" - MIT Press
% calculate the size of the wavelet
n_wavelet = length(wavelet);
half_of_wavelet_size = round((n_wavelet-1) / 2);
% number of data points in the input data
numTimePoints = length(inputData);
% together these determine the number of wavelet convolutions
n_convolution = n_wavelet+numTimePoints-1;

% Run Fourier transform
fft_EEG = fft(inputData,n_convolution);

% Fourier transform of the wavelet for calculation
fft_wavelet = fft(wavelet,n_convolution);
% convolve wavelet with data and inverse Fourier transform back out
convolution_result_fft = ifft(fft_wavelet.*fft_EEG,n_convolution);
% Pull the actual data from the center of this convolution
if mod(n_wavelet,2)==1
    % odd length wavelet
    evenFlag = 0;
else
    evenFlag = 1;
end
convolution_result_fft = convolution_result_fft((half_of_wavelet_size+1-evenFlag):(end-half_of_wavelet_size));

% Calculate the phase of the resulting convolution
freqPhase = angle(convolution_result_fft);
% Calculate the amplitude of the resulting convolution
freqAmp = abs(convolution_result_fft);
                        
end % end of function