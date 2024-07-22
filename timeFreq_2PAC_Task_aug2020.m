% Updated version for August 14th, 2020
%   Runs logarithmically spaced oscillations


% Goal of the script - run single subject analysis
%   1. Run time frequency analysis

clc; clear all; close all;

%% Variables that are liable to change
preprocVersion = 'oct2018';
timeFreqVersion = 'aug2020';

chan128 = 1;
if chan128
    chan128_str = '_128';
else
    chan128_str = '';
end


%% Paths to Toolboxes
% Root is specific to the computer you are working on
COMPUTER = '/Users/justinri/';
ROOT = [COMPUTER 'Dropbox (Frohlich Lab)/'];
CODE = [ROOT 'Frohlich Lab Team Folder/Codebase/CodeJustin/'];

% Riddler Toolbox
RIDDLER_TOOLBOX = [CODE 'Riddler_Toolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));

% Toolboxes
TOOLBOXES = [CODE 'Toolboxes/'];
EEGLAB_TOOLBOX = [TOOLBOXES 'eeglab14_0_0b/'];
addpath(EEGLAB_TOOLBOX);
[ALLEEG,EEG,CURRENTSET] = eeglab('rebuild');

%% Experiment paths
% All data
DATA = [ROOT 'HumanStudies/2PAC/'];
% Preprocessed EEG data
PROC_EEG = [DATA 'ProcEEG/'];
% Individual level of analysis
INDIV_EEG = [DATA 'Individual_AnalysisEEG/'];
mkdir_JR(INDIV_EEG);
% Time Frequency directory for individual subject data
INDIV_TF = [INDIV_EEG 'TimeFreq_Task_' timeFreqVersion '/'];
mkdir_JR(INDIV_TF);
% Individual level of analysis
GROUP_EEG = [DATA 'GroupEEG/'];
mkdir_JR(GROUP_EEG);
% Time Frequency directory for individual subject data
GROUP_TF = [GROUP_EEG 'TimeFreq_Task_' timeFreqVersion '/'];
mkdir_JR(GROUP_TF);

%% Subjects
SUBJECTS = {...
    '01','02','03','04',...
    '05','06','07','08',...
    '09','10','11','12',...
    '14','15','16',... % missing sub 13 - non-compliant in task
    '17','18','19','20',...
    '21','22','23','24'};
numSub = length(SUBJECTS);

%% Task info
COND = {'R4','R8','D1','D2'};
numCond = length(COND);

%% Load EGI channel location information
chanlocFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_channelLocations.mat'];
chanlocsStruct = load(chanlocFile);
chanlocs = chanlocsStruct.chanlocs;
dataChannelsFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_dataChannels_idxs.mat'];
dataChannelsStruct = load(dataChannelsFile);
dataChannels = dataChannelsStruct.dataChannels;
numDataChannels = length(dataChannels);
data_chanlocs = chanlocs(dataChannels);

%% Time Frequency Information
% Relevant PAC info
% Frequency range specifications
lowerBoundFreq = 2;
upperBoundFreq = 58;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
[freq] = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);
numFreqs = length(freq);
numCycles = 5;
sampleRate = 200;
preStim = -1;
postStim = 2;
analysisWindow = postStim - preStim;
ms = 1000;
timePoints = (preStim*ms):(ms/sampleRate):((postStim*ms)-(ms/sampleRate));
numTimePoints = length(timePoints);
%mirrorTimePoints = timePoints + (analysisWindow *ms);

%% Single Subject Time Frequency analysis
timeFreqFiles = cell(1,numSub);
% Loop through each subject
for subIdx = 1:numSub
    subject = ['sub' SUBJECTS{subIdx}];
    day = 'Baseline';
        
    % Frequently used for naming sysstem
    daySub = [subject '_' day];
    % Prefix for use in calculating TF
    timeFreqFileprefix = ['tf_' daySub '_' timeFreqVersion];
    % output file for time frequency results
    timeFreqFile = sprintf('%s%s_%03dchannels_%ifreqs_%icycles.mat',...
        INDIV_TF,timeFreqFileprefix,numDataChannels,numFreqs,numCycles);
    timeFreqFiles{subIdx} = timeFreqFile;
    % Check if time freq file already exists
    if exist(timeFreqFile,'file')~=2
        fprintf('\n\nProcessing time frequency for %s %s\n',subject, day);
        % Gather preprocessed EEG data
        SUB_PROC = [PROC_EEG subject '/Baseline_' preprocVersion '/'];
        % Initialize output variables
        allCond_ERSP = NaN(numFreqs,numTimePoints,numDataChannels,numCond);
        allCond_ITPC = NaN(numFreqs,numTimePoints,numDataChannels,numCond);
        numCondExist = 0;
        
        COND_EEG_FILES = cell(1,numCond);
        % Loop through each in condition
        for condIdx = 1:numCond
            cond = COND{condIdx};
            
            % Load preprocessed epoched data
            preprocFile = [SUB_PROC 'final_ica_' cond '_labeled_preproc_2PAC_' daySub '.set'];
            assert(exist(preprocFile,'file')==2);
            
            % Load previous data and mirror flip
            COND_EEG_FILES{condIdx} = preprocFile;
        end
        
        tfInfo = struct(...
            'conditions',{COND},...
            'sampleRate',sampleRate,...
            'numCycles',numCycles,...
            'fileName',timeFreqFileprefix,...
            'frequencies',freq,...
            'mirrorFlag',1,...
            'dataChannels',dataChannels,...
            'baselinePeriod',[0.4 0.7]);
            
        % Run time frequency analysis
        [tfOutputFiles,~] = timeFreq_JR(INDIV_TF,COND_EEG_FILES,tfInfo);
    end % time frequency files exists
end % loop subjects

%% Group level time frequency analysis
ANALYSES = {...
    'Abstraction','D1_D2_R4_R8';...
    'SetSize',    'R8_D2_R4_D1';...
    'SingleR4',   'R4_0'  ;...
    'SingleR8',   'R8_0'  ;...
    'SingleD1',   'D1_0'  ;...
    'SingleD2',   'D2_0'  ;...
    'D2vsR8',     'D2_R8' };
numAnalyses = size(ANALYSES,1);
numColumns = 2;

% Run permutation testing in subset of electrodes
ROIS = {...
    'FCz',6;...
    'aPFC',[ 11   4   5  10  12  16  18  19]};
ROI_NAMES = ROIS(:,1);
ROI_IDXS = ROIS(:,2);
numROIs = length(ROI_NAMES);
numIterations = 1000;

% Two types of analysis
ANALYSIS_TYPE = [{'eachChannel'} ROI_NAMES'];
% Each channel is going to just calculate the contrast: t,p,db for each
% channel. The second ROI analysis will run permutation testing for the
% electrodes selected

% Apply a threshold to the cluster size for calculating mass
% otherwise spurious tiny clusters will warp this
minClusterSize = 30;
alphaThreshold = 0.05;

% Loop analyses
for analysisIdx = 1:numAnalyses
    analysisName = ANALYSES{analysisIdx,1};
    
    clear analysisData
    
    % All the data needed to calculate output above
    analysisData = NaN(numTimePoints,numFreqs,numSub,numColumns);

    % Analysis data is generated or loaded on the first pass
    for typeIdx = 1:length(ANALYSIS_TYPE)
        analysisType = ANALYSIS_TYPE{typeIdx};
        if strcmp(analysisType,'eachChannel')
            roiFLAG = 0;
            channelIdxs = 1:numDataChannels;
            numDataOut = numDataChannels;
        else
            roiFLAG = 1;
            roiIdx = typeIdx - 1;
            channelIdxs_128 = ROI_IDXS{roiIdx};
            [~,channelIdxs,~]=intersect(dataChannels,channelIdxs_128);
            numDataOut = 1;
            analysisType = sprintf('%s_alpha%03d_size%i',analysisType,alphaThreshold*1000,minClusterSize); 
        end
        
        % Analysis output file
        analysisOutputFile = [GROUP_TF analysisName '_' analysisType '_' timeFreqVersion '.mat'];
        if exist(analysisOutputFile,'file')~=2
            
            fprintf('\n\nProcessing %s %s\n',analysisName,analysisType);

            % Analysis logic
            analysisLogic = regexp(ANALYSES{analysisIdx,2},'_','split');
            numLogic = length(analysisLogic);
            halfLogic = numLogic / 2;
            columns = {...
                analysisLogic(1:halfLogic),...
                analysisLogic((halfLogic+1):end)};

            % Gather data before averaging
            colData = zeros(numTimePoints,numFreqs,numDataOut,numSub,numColumns);
            
            % Loop through two columns
            for columnIdx = 1:numColumns
                column = columns{columnIdx};
                numCond_col = length(column);
                condIdxs = NaN(1,numCond_col);
                for condColIdx = 1:numCond_col
                    foundCondIdx = find(strcmpi(column{condColIdx},COND));
                    if ~isempty(foundCondIdx)
                        condIdxs(condColIdx) = foundCondIdx;
                    end
                end
                % skip this column if no conditions - for single "0" column
                if ~isempty(find(isnan(condIdxs),1))
                    continue;
                end
                
                % Loop subjects
                for subIdx = 1:numSub
                    timeFreqFile = timeFreqFiles{subIdx};
                    timeFreq = load(timeFreqFile);
                    % ERSP: time x freq x chan x cond
                    
                    erspCol = NaN(numTimePoints,numFreqs,numDataOut,numCond_col);
                    for condColIdx = 1:numCond_col
                        condIdx = condIdxs(condColIdx);
                        cond = COND{condIdx};
                        condTF = timeFreq.tf.(cond).ersp;
                        chanCondTF = condTF(:,:,channelIdxs);
                        if length(channelIdxs)>1
                            chanCondTF = mean(chanCondTF,3);
                        end
                        erspCol(:,:,:,condColIdx) = chanCondTF;
                    end
                    
                    % Average across conditions for this column
                    if numCond_col > 1
                        erspCol = squeeze(mean(erspCol,4));
                    end
                    % If doing an ROI permutation analysis, then average
                    % across electrodes in that ROI
                    if roiFLAG && length(channelIdxs)>1
                        erspCol = squeeze(mean(erspCol,3));
                    end
                    % Store data for this column
                    colData(:,:,:,subIdx,columnIdx) = erspCol;
                end % loop subjects
            end % loop two columns
            
            
            % Output is a t-test & db
            outDim = [numTimePoints numFreqs numDataOut];
            analysisOut_t = NaN(outDim);
            analysisOut_p = NaN(outDim);
            analysisOut_db = NaN(outDim);
            
            % Loop through data and run analysis
            for timeIdx = 1:numTimePoints
                for freqIdx = 1:numFreqs
                    for dataIdx = 1:numDataOut
                        
                        % Compare two conditions
                        col1 = squeeze(colData(timeIdx,freqIdx,dataIdx,:,1));
                        col2 = squeeze(colData(timeIdx,freqIdx,dataIdx,:,2));
                        [~,p,~,tstat]=ttest(col1,col2);
                        
                        % Sometimes column 2 is zeros
                        colDiff = mean(col1 - col2);
                        
                        % Store analysis out
                        analysisOut_t(timeIdx,freqIdx,dataIdx) = tstat.tstat;
                        analysisOut_p(timeIdx,freqIdx,dataIdx) = p;
                        analysisOut_db(timeIdx,freqIdx,dataIdx) = colDiff;
                    end % loop data
                end % loop freq
            end % loop time
            
            
            % Save out tstat, p-values, and decibels for all data channels
            analysisStruct_real = struct(...
                't',analysisOut_t,...
                'p',analysisOut_p,...
                'db',analysisOut_db);
            
            % Analysis Structure to save out
            analysisStruct = struct(...
                'real',analysisStruct_real);
            
            if roiFLAG
                % Run permutation testing on data
                iter_maxSize = NaN(1,numIterations);
                iter_maxT = NaN(1,numIterations);
                iter_maxMass = NaN(1,numIterations);
                
                % No channels for ROI analysis - just one ROI
                colData = squeeze(colData);
                
                % Loop through iterations
                fprintf('iter ');
                for iterIdx = 1:numIterations
                    
                    % update the user
                    iterStr = fprintf('%04d',iterIdx);
                    
                    % Randomize subject assignment
                    randSub = randperm(numSub);
                    if mod(numSub,2) ~= 0
                        extraSub = randi(2) - 1;
                    else
                        extraSub = 0;
                    end
                    halfSub = floor(numSub/2) + extraSub;
                    randSwapSub = randSub(1:halfSub);
                    
                    % Gather the t-stat for this iteration
                    iterData_t = NaN(numTimePoints,numFreqs);
                    iterData_p = NaN(numTimePoints,numFreqs);
                    
                    % Loop through data and run analysis
                    for timeIdx = 1:numTimePoints
                        for freqIdx = 1:numFreqs
                            % Compare two conditions
                            col1 = colData(timeIdx,freqIdx,:,1);
                            col2 = colData(timeIdx,freqIdx,:,2);

                            % randomly swap half of the subjects
                            tmp = col1;
                            col1(randSwapSub) = col2(randSwapSub);
                            col2(randSwapSub) = tmp(randSwapSub);

                            % run the t-test on random data
                            [~,p,~,tstat]=ttest(col1,col2);

                            % Store analysis out
                            iterData_t(timeIdx,freqIdx) = tstat.tstat;
                            iterData_p(timeIdx,freqIdx) = p;
                            
                        end % loop freq
                    end % loop time
                    
                    [maxSize,maxT,maxMass] = significanceTimeFreq_JR(...
                        iterData_t,iterData_p,alphaThreshold,minClusterSize,'test');
                    
                    % Calculate max size, t, and mass
                    iter_maxSize(iterIdx) = maxSize;
                    iter_maxT(iterIdx)    = maxT;
                    iter_maxMass(iterIdx) = maxMass;
                    
                    % update the user
                    fprintf('\b\b\b\b');
                end % loop iterations
                
                % Significance index with alpha = 0.05
                sigIdx = ceil(numIterations * alphaThreshold);
                
                % Sort iterations
                iter_maxSize = sort(iter_maxSize,'descend');
                iter_maxT = sort(iter_maxT,'descend');
                iter_maxMass = sort(iter_maxMass,'descend');
                
                % Get significance threshold
                sig_size = iter_maxSize(sigIdx);
                sig_maxT = iter_maxT(sigIdx);
                sig_mass = iter_maxMass(sigIdx);
                
                % Store output from the permutation testing
                sigStruct = struct(...
                    'alpha',alphaThreshold,...
                    'minSize',minClusterSize,...
                    'size',sig_size,...
                    't',sig_maxT,...
                    'mass',sig_mass);
                iterStruct = struct(...
                    'numIterations',numIterations,...
                    'size',iter_maxSize,...
                    't',iter_maxT,...
                    'mass',iter_maxMass);
                permutationStruct = struct(...
                    'sig',sigStruct,...
                    'iterations',iterStruct);
                
                % Gather relevant information
                analysisStruct.permutation = permutationStruct;
            end % if running permutation analysis
            
            % Save analysis output
            save(analysisOutputFile,'-struct','analysisStruct');
            
        end % loop analysis type
    end % analysis file exists
end % loop analyses