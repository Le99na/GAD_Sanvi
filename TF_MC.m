% wipe out all previous data
clear all; clc; close all;

% Preprocessing version
preprocVersion = '02052023';

% Use ICLabel or not
useIClabel = 1;
skipManualIC = 0;

%% Step 0: get organized
% Paths to data
root = 'C:\Users\magcam\Documents\code_matlab\Time_frequency_analysis\'
addpath(genpath(root))
ClassObj = NNC_path();

root_raw = ClassObj.DataRoot_GAD_anxious_original
root_proc = ClassObj.DataRoot_GAD_anxious_processed

%generate cell array of available participants
grps_strct=dir(root_raw)
grp_nams = {grps_strct.name}

longerThan3 = cellfun(@(x) length(x) > 3 && contains(x,'ptp'), grp_nams)
GROUPS = grp_nams(find(longerThan3))

numGroups = length(GROUPS);

% Add the EEGLAB toolbox to the path
addpath(genpath(ClassObj.EEGLAB_TOOLBOX));

%add csc toolbox
addpath(genpath(ClassObj.CSC_TOOLBOX));

%add epi-eeg toolbox
addpath(genpath(ClassObj.epi_TOOLBOX));

%add power scripts toolbox
addpath(genpath(ClassObj.pwr_TOOLBOX));

% Add Mengsen butter function + rainbow maps

addpath(genpath(ClassObj.rainbow1)); 

addpath(genpath(ClassObj.rainbow2));

% FREQ = [1:100]
% numFreqs = length(FREQ);

sampleRate = 200;
numCycles = 3; %high temporal resolution: 3-5, high frequency resolution: 5-8
baselinePeriod = 1:40
lowerBoundFreq = 2;
upperBoundFreq = 40;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
FREQ = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);

% Number of cycles to estimate tf can change for each frequency
if length(numCycles) ~= numFreqs
    numCycles = ones(1,numFreqs) * numCycles(1);
end

% Defualt is to mirror the data
mirrorData_FLAG = 1;


% Start EEGlab
eeglab

%% Go through each participant
for groupIdx = 1:length(fileList)
    % Name of the group which is used in folder organization
    % path_in1 = [root_raw GROUPS{groupIdx} '\EEG\']
    % fls = {dir(path_in1).name}
    % longerThan3 = fls{cellfun(@(x) length(x) > 3,fls)}
    % path_in2 = [path_in1 longerThan3 '\']

    path_out = [root_proc GROUPS{groupIdx} '_nc3_params\']

    if exist(path_out,'dir')~=7
        mkdir(path_out);
    end
    
    %read-in relevant file
    % fls_all = {dir(path_in2).name}
    % fl_rel_idx = find(contains(fls_all,'ERN_vals.set'))
    % fl_rel = fls_all{fl_rel_idx}

    EEG = pop_loadset(fileList{groupIdx})

    answ = {EEG.event.type}
    idx_ep_cor = []
    idx_ep_icor = []
    count=1 %this is the count for the number of epochs
    for i = 1:numel(answ)
        if all(contains(answ{i},'PR_correct'))
            i
            count
            idx_ep_cor = [idx_ep_cor,count];
            count=count+1;
            answ{i}
        elseif all(contains(answ{i},'PR_incorrect'))
            i
            count
            idx_ep_icor = [idx_ep_icor,count];
            count=count+1;
            answ{i}
        end
    end

    EEG_corr = EEG.data([6,11,129],:,idx_ep_cor); %6:FCz, 11:Fz, 129:Cz
    EEG_icorr = EEG.data([6,11,129],:,idx_ep_icor);

    EEG.idx_cor = idx_ep_cor
    EEG.idx_icor = idx_ep_icor

    EEG = pop_saveset(EEG, 'filename', [GROUPS{groupIdx}, '_idx_cor_icor.set'], 'filepath', path_out);

    roi_Data_cor = squeeze(mean(EEG_corr(:,:,:),1));

    roi_Data_icor = squeeze(mean(EEG_icorr(:,:,:),1));

    % find the save of this condition's data
    [numTimePts_cor,numTrials_cor] = size(roi_Data_cor);
    [numTimePts_icor,numTrials_icor] = size(roi_Data_icor);

    % Run ersp for each channel separately
    amp_cor = NaN(numTimePts_cor,numFreqs,numTrials_cor);
    amp_icor = NaN(numTimePts_icor,numFreqs,numTrials_icor);
    phase_cor = NaN(numTimePts_cor,numFreqs,numTrials_cor);
    phase_icor = NaN(numTimePts_icor,numFreqs,numTrials_icor);
  
    % loop through each trial for this subject
    for trialIdx = 1:numTrials_cor
        % data for only this trial
        trialData = roi_Data_cor(:,trialIdx)';
        
        % If the data is being mirror flipped
        if mirrorData_FLAG
            % mirror flip the trial data
            mirrorData = [fliplr(trialData) trialData fliplr(trialData)];
            mirrorIdxs = (numTimePts_cor+1):(numTimePts_cor*2);
        else
            % Or just use the original data
            mirrorData = trialData;
            mirrorIdxs = 1:numTimePts_cor;
        end
        
        % run time frequency analysis for each frequency
        for freqIdx = 1:numFreqs
            freq = FREQ(freqIdx);
            % run wavelet based time frequency estimate
            [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles(freqIdx));
            
            % Prepare output
            phase_cor(:,freqIdx,trialIdx) = phase(mirrorIdxs);
            amp_cor(:,freqIdx,trialIdx) = amp(mirrorIdxs);
        end % freq loop
    end % trial loop

    % loop through each trial for this subject
    for trialIdx = 1:numTrials_icor
        % data for only this trial
        trialData = roi_Data_icor(:,trialIdx)';
        
        % If the data is being mirror flipped
        if mirrorData_FLAG
            % mirror flip the trial data
            mirrorData = [fliplr(trialData) trialData fliplr(trialData)];
            mirrorIdxs = (numTimePts_icor+1):(numTimePts_icor*2);
        else
            % Or just use the original data
            mirrorData = trialData;
            mirrorIdxs = 1:numTimePts_icor;
        end
        
        % run time frequency analysis for each frequency
        for freqIdx = 1:numFreqs
            freq = FREQ(freqIdx);
            % run wavelet based time frequency estimate
            [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles(freqIdx));
            
            % Prepare output
            phase_icor(:,freqIdx,trialIdx) = phase(mirrorIdxs);
            amp_icor(:,freqIdx,trialIdx) = amp(mirrorIdxs);
        end % freq loop
    end % trial loop

    % Average across trials
    amp_trialAvg_cor = mean(amp_cor,3);
    amp_trialAvg_icor = mean(amp_icor,3);
    
    % Baseline correction
    amp_baselineCorrected_cor = NaN(numTimePts_cor,numFreqs);
    amp_baselineCorrected_icor = NaN(numTimePts_icor,numFreqs);

    % Calculate baseline
    amp_baselineAvg_cor = mean(amp_trialAvg_cor(baselinePeriod,:),1);
    amp_baselineAvg_icor = mean(amp_trialAvg_icor(baselinePeriod,:),1);

    % Baseline correct all data
    for timeIdx = 1:numTimePts_cor
        timepoint_amp = amp_trialAvg_cor(timeIdx,:);
        timepoint_amp = (timepoint_amp - amp_baselineAvg_cor) ./ amp_baselineAvg_cor;
        amp_baselineCorrected_cor(timeIdx,:) = timepoint_amp;

        timepoint_amp = amp_trialAvg_icor(timeIdx,:);
        timepoint_amp = (timepoint_amp - amp_baselineAvg_icor) ./ amp_baselineAvg_icor;
        amp_baselineCorrected_icor(timeIdx,:) = timepoint_amp;

    end

    EEG.amp_avg_cor = amp_baselineCorrected_cor
    EEG.amp_avg_icor = amp_baselineCorrected_icor

    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_avg_amp.set'], 'filepath', path_out);

    % Plot the spectrogram
    figure;
    imagesc(1/200:1/200:300/200, FREQ, amp_baselineCorrected_cor.');
    set(gca, 'YDir', 'normal');
    title([EEG.filename(1:6) ' time-frequency spectrogram cor']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-5 5])
    colormap(flipud(brewermap([],'RdBu')));
    colorbar;
    % Save the figure as a JPEG
    output_filename = [path_out EEG.filename(1:6) '_time_frequency_spectrogram_cor.jpg'];  % Specify the output filename
    saveas(gcf, output_filename, 'jpeg');
    close all

    % Plot the spectrogram
    figure;
    imagesc(1/200:1/200:300/200, FREQ, amp_baselineCorrected_icor.');
    set(gca, 'YDir', 'normal');
    title([EEG.filename(1:6) ' time-frequency spectrogram icor']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    % Save the figure as a JPEG
    output_filename = [path_out EEG.filename(1:6) '_time_frequency_spectrogram_icor.jpg'];  % Specify the output filename
    saveas(gcf, output_filename, 'jpeg');
    close all
    
end



