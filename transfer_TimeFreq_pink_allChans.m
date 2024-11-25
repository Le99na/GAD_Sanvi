% wipe out all previous data
clear all; clc; close all;

scripts = '/work/users/m/a/magcam/GAD/scripts/';
addpath(scripts)

root_raw = '/work/users/m/a/magcam/GAD/EEG_mats_all/';

% Define the top-level directory and the pattern
topLevelDir = root_raw; % Replace with your directory path
filePattern = '**/EEG*.mat'; % Replace with your file pattern (e.g., '*.txt' or '*.mat')

% Search for all matching files
files = dir(fullfile(topLevelDir, filePattern));
% Preallocate cell array to store full paths
filePaths = cell(length(files), 1);
% Store the full paths of the matching files in the cell array
for k = 1:length(files)
    filePaths{k} = fullfile(files(k).folder, files(k).name);
end

sampleRate = 200;
numCycles = 5

lowerBoundFreq = 2;
upperBoundFreq = 58;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
FREQ = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);

numFreq = length(FREQ);
baselinePeriod = 1:40;

% Defualt is to mirror the data
mirrorData_FLAG = 1;

numGroups = length(filePaths);


for subIdx = 1:numGroups
    fl = filePaths{subIdx};
    load(fl);
    dat_c = EEG.data(:,:,EEG.idx_corr);
    dat_i = EEG.data(:,:,EEG.idx_icorr);

    prts = split(fl,'/');
    prts2 = split(prts{end},'_');
    id = prts2{end}(4:6); 

    pwr_i = zeros(129,300,size(dat_i,3),numFreq);
    phs_1 = zeros(129,300,size(dat_i,3),numFreq);

    for chn = 1:129

        for trlIdx = 1:size(dat_i,3)
            i_Data = dat_i(chn,:,trlIdx);
    
            % find the size of this condition's data
            numTimePts = size(i_Data,2);
    
            amp_i = NaN(length(i_Data),numFreq);
            phase_i = NaN(length(i_Data),numFreq);
    
            % If the data is being mirror flipped
            if mirrorData_FLAG
                % mirror flip the trial data
                mirrorData = [fliplr(i_Data) i_Data fliplr(i_Data)];
                mirrorIdxs = (numTimePts+1):(numTimePts*2);
            else
                % Or just use the original data
                mirrorData = dat_i;
                mirrorIdxs = 1:length(dat_i);
            end
    
            for freqIdx = 1:numFreq
                freq = FREQ(freqIdx);
                % run wavelet based time frequency estimate
                [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles);
                
                % Prepare output
                phase_i(:,freqIdx) = phase(mirrorIdxs);
                amp_i(:,freqIdx) = amp(mirrorIdxs);
            end
    
            % Baseline correction
            amp_baselineCorrected_i = NaN(numTimePts,numFreq);
        
            % Calculate baseline
            amp_baselineAvg_i = mean(amp_i(baselinePeriod,:),1);
        
            % Baseline correct all data
            for timeIdx = 1:numTimePts
                timepoint_amp = amp_i(timeIdx,:);
                timepoint_amp = (timepoint_amp - amp_baselineAvg_i) ./ amp_baselineAvg_i;
                amp_baselineCorrected_i(timeIdx,:) = timepoint_amp;
            end
    
            pwr_i(chn,:,trlIdx,:) = amp_baselineCorrected_i;
            phs_i(chn,:,trlIdx,:) = phase_i;
    
        end
    end

    save(['/work/users/m/a/magcam/GAD/time_freq_perTrial_andChan_pink/pwr_Sub' id '_incorrect' num2str(trlIdx) 'pink.mat'], 'pwr_i','-v7.3');
    save(['/work/users/m/a/magcam/GAD/time_freq_perTrial_andChan_pink/phs_Sub' id '_incorrect' num2str(trlIdx) 'pink.mat'], 'phs_i','-v7.3');


    pwr_c = zeros(129,300,size(dat_c,3),numFreq);
    phs_c = zeros(129,300,size(dat_c,3),numFreq);

    for chn = 1:129

        for trlIdx = 1:size(dat_c,3)
            c_Data = dat_c(chn,:,trlIdx);
    
            % find the size of this condition's data
            numTimePts = size(c_Data,2);
    
            amp_c = NaN(length(c_Data),numFreq);
            phase_c = NaN(length(c_Data),numFreq);
    
            % If the data is being mirror flipped
            if mirrorData_FLAG
                % mirror flip the trial data
                mirrorData = [fliplr(c_Data) c_Data fliplr(c_Data)];
                mirrorIdxs = (numTimePts+1):(numTimePts*2);
            else
                % Or just use the original data
                mirrorData = dat_c;
                mirrorIdxs = 1:length(dat_c);
            end
    
            for freqIdx = 1:numFreq
                freq = FREQ(freqIdx);
    
                % run wavelet based time frequency estimate
                [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles);
                
                % Prepare output
                phase_c(:,freqIdx) = phase(mirrorIdxs);
                amp_c(:,freqIdx) = amp(mirrorIdxs);
            end
    
            % Baseline correction
            amp_baselineCorrected_c = NaN(numTimePts,numFreq);
        
            % Calculate baseline
            amp_baselineAvg_c = mean(amp_c(baselinePeriod,:),1);
        
            % Baseline correct all data
            for timeIdx = 1:numTimePts
                timepoint_amp = amp_c(timeIdx,:);
                timepoint_amp = (timepoint_amp - amp_baselineAvg_c) ./ amp_baselineAvg_c;
                amp_baselineCorrected_c(timeIdx,:) = timepoint_amp;
            end
    
            pwr_c(chn,:,trlIdx,:) = amp_baselineCorrected_c;
            phs_c(chn,:,trlIdx,:) = phase_c;
    
        end
    end

    save(['/work/users/m/a/magcam/GAD/time_freq_perTrial_andChan_pink/pwr_Sub' id '_correct_' num2str(trlIdx) 'pink.mat'], 'pwr_c','-v7.3');
    save(['/work/users/m/a/magcam/GAD/time_freq_perTrial_andChan_pink/phs_Sub' id '_correct_' num2str(trlIdx) 'pink.mat'], 'phs_c','-v7.3');

end

