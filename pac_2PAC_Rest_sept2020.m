function pac_2PAC_Rest_sept2020(subIdx,varargin)

% Important options
highDefinitionFLAG = 0;

% default is to be ready for the longleaf cluster
if ~isempty(varargin)
    localComputer = varargin{1};
else
    localComputer = 0;
end

analysisVersion = 'aug2020';

% Subjects
SUBJECTS = {...
    '01','02','03','04',...
    '05','06','07','08',...
    '09','10','11','12',...
    '14','15','16',... '13' non-compliant subject
    '17','18','19','20',...
    '21','22','23','24'};
numSub = length(SUBJECTS);

subject = ['sub' SUBJECTS{subIdx}];

% Version of data
preprocVersion = 'apr2019';

if localComputer
    %% LOCAL COMPUTER - RUNNING ON DROPBOX
    % root
    %ROOT = 'C:/Users/Justin/Dropbox (Frohlich Lab)/';
    ROOT = '/Users/justinri/Dropbox (Frohlich Lab)/';

    % Toolboxes
    CODE = [ROOT 'Frohlich Lab Team Folder/Codebase/CodeJustin/'];
    
    % Only load this locally
    EEGLAB_TOOLBOX = [CODE 'Toolboxes/eeglab14_1_2b/'];
    addpath(EEGLAB_TOOLBOX);
    eeglab('rebuild');
    
    % Data paths
    DATA = [ROOT 'HumanStudies/2PAC/'];
    RAW_BEHAV = [DATA 'RawBehavior/'];
    PROC_EEG = [DATA 'ProcEEG/'];
    GROUP_EEG = [DATA 'GroupEEG/'];
    mkdir_JR(GROUP_EEG);
    
    GROUP_PAC = [GROUP_EEG 'PAC_Rest_' analysisVersion '/'];
    mkdir_JR(GROUP_PAC);
    SINGLE_PAC = sprintf('%sIndividual_AnalysisEEG/PAC_Rest_%s/',DATA,analysisVersion);
    mkdir_JR(SINGLE_PAC);

    %% Loop through participants and calculate PAC
    subject =['sub' SUBJECTS{subIdx}];

    % output file for single subject level
    SUB_PAC_OUTPUT = [SINGLE_PAC subject '/'];
    mkdir_JR(SUB_PAC_OUTPUT);
    
    SUB_PROC_EEG  = [PROC_EEG subject '/'];
    
    SUB_RAW_BEHAV = [RAW_BEHAV subject '/'];
    
else
    
    %% LONGLEAF CLUSTER
    % NAS Root
    NAS_ROOT = '/nas/longleaf/home/justinri/';
    % Pine Root
    PINE_ROOT = '/pine/scr/j/u/justinri/2PAC/';

    % Code base on cluster
    CODE = [NAS_ROOT 'CodeJustin/'];
    
    % Data paths
    PROC_EEG_IN = [PINE_ROOT 'ProcEEG_Rest/'];
    
    BEHAV_IN = [PINE_ROOT 'RawBehavior/'];

    % output file for single subject level
    PAC_OUTPUT =[PINE_ROOT 'PAC_Rest/'];
     
end

% Load riddler toolbox no matter what
RIDDLER_TOOLBOX = [CODE 'Riddler_Toolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));

% Task information
COND = {'R4','R8','D1','D2'};
numCond = length(COND);

% Contrasts of interest
CONTRASTS = {...
    'Abstraction',...
    'Setsize',...
    'SingleR4',...
    'SingleR8',...
    'SingleD1',...
    'SingleD2'};
% Logic: split in half and add each half, then compare them
CONTRAST_LOGIC = {...
    [3 4],[1 2];...
    [2 4],[1 3];...
    1,0;...
    2,0;...
    3,0;...
    4,0};
numContrasts = length(CONTRASTS);

% Task timing information
sampleRate = 200;
numTimePoints = 600;
numSeconds = 3.0;
sampleWindow = numSeconds / numTimePoints;
realTime = -1:sampleWindow:(2-sampleWindow);

%% Load EGI channel location information
dataChannelsFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_dataChannels_idxs.mat'];
dataChannelsStruct = load(dataChannelsFile);
dataChannels = dataChannelsStruct.dataChannels;
numDataChannels = length(dataChannels);

% Region of interest look up table
ROI_DEF = {...
    'FCz',   6                                ,... % centerFCz (e6)
    'aPFC',  [ 11   4   5  10  12  16  18  19],... % center Fz (e11)
    'lPFC',  [ 19  20  23  24  27  28]        ,... % center F3 (e24)
    'rPFC',  [  3   4 117 118 123 124]        ,... % center F4 (e124)
    'cPFC',  [  5   6   7  12  13 106 112]    ,... % center FCz (e6)
    'lM1',   [ 36  29  30  35  37  41  42]    ,... % center C3 (e36)
    'rM1',   [104  87  93 103 105 110 111]    ,... % center C4 (e104
    'lPar',  [ 52  42  47  51  53  59  60]    ,... % center P3 (e52)
    'rPar',  [ 92  85  86  91  93  97  98]    ,... % center P4 (e92)
    'ParOcc',[ 61  62  67  72  77  78]        ,... % center POz (e62)
    'Occ',   [ 75  70  71  76  83]            ,... % center Oz (e75)
    'rParOcc',[62 72 76 77 78 83 84 85 90 91 92],...
    'aParOcc',[37 42 53 54 61 78 79 86 87 93]}; % derived from theta-gamma coupling
numROIs = length(ROI_DEF) / 2;
ROI_NAMES = ROI_DEF(1:2:length(ROI_DEF));

%% Version 1 - Hilbert transform in canonical bands
delta  = [ 2  3];
theta  = [ 4  7];
beta   = [15 25];
lowBeta= [10 20];
gamma  = [35 58];
% Two specific hypotheses
PAC_PAIRS = {'delta_beta','theta_gamma'};
numPacPairs = length(PAC_PAIRS);
canonical_SEED_ROIS = {};%'cPFC','aPFC'};%,'rPFC'};
numSeedROIs_canonical = length(canonical_SEED_ROIS);

if highDefinitionFLAG
    
    % 1/f scaled time frequency estimate
    % high definition
    numFreqs = 150;

else
    % 1/f scaled time frequency estimate
    % low definition
    numFreqs = 75;
end

%% Exhaustive PAC
% Relevant PAC info
% Frequency range specifications
lowerBoundFreq = 2;
upperBoundFreq = 58;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
[freq] = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);
% Freq range for PAC
lowFreqRange  = [2 8];
highFreqRange = [9 58];
lowFreqIdx_start  = find(lowFreqRange(1)  < (freq+0.001),1);
lowFreqIdx_stop   = find(lowFreqRange(2)  < (freq+0.001),1);
highFreqIdx_start = find(highFreqRange(1) < (freq+0.001),1);
highFreqIdx_stop  = find(highFreqRange(2) < (freq+0.001),1);
lowFreqIdxs  = lowFreqIdx_start:lowFreqIdx_stop;
highFreqIdxs = highFreqIdx_start:highFreqIdx_stop;
LOW_FREQ  = freq(lowFreqIdxs);
HIGH_FREQ = freq(highFreqIdxs);
numLowFreq  = length(LOW_FREQ);
numHighFreq = length(HIGH_FREQ);
numCycles_low = 3;
numCycles_high = 5;
exhaustive_SEED_TARGET = {'aPFC_aParOcc','cPFC_rM1'};%'cPFC_cPFC','aPFC_aPFC'};
numPairsROI_exhaustive = length(exhaustive_SEED_TARGET);

% Iterations for PAC analysis
numIterations = 1000;
% Number of phase bins
numPhaseBins = 30;
pacOptions = struct(...
    'numIterations',numIterations,...
    'numPhaseBins',numPhaseBins);

% Method for calculating PAC
METHODS = {'cohen','tort','canolty',};
numMethods = length(METHODS);

% Final version of preprocessing
preprocPrefix = 'tmcgabrdfei_2PAC_';

%% Stimulation Day information
STIM_DAYS = {'Baseline','Stim1','Stim2','Stim3'};
numDays = length(STIM_DAYS);
TACS_CONDITIONS = {'Baseline','deltaBeta','thetaGamma','sham'};
numTacsConditions = length(TACS_CONDITIONS);

% Information on task blocks
numBlocks = 8;
numRestPeriods = numBlocks + 1;

% Loop through days
for dayIdx = 1:numDays
    day = STIM_DAYS{dayIdx};

    if localComputer
        % Subject and day raw behavior
        BEHAV_IN = [SUB_RAW_BEHAV day '/'];
        
        % Preprocessing directory
        PROC_EEG_IN = [SUB_PROC_EEG day '_Rest_' preprocVersion '/'];
        
        % PAC output directory
        PAC_OUTPUT = SUB_PAC_OUTPUT;
        mkdir_JR(PAC_OUTPUT);
        
    end

    % Unblind to the condition of stimulation
    if strcmp(day,'Baseline')
        stimCondition = 'Baseline';
    else
        % Each block has a results file
        taskBlock_resultsFileFndr = dir([BEHAV_IN 'results_2PAC_' subject '_' day '_1_*.mat']);
        taskBlock_resultsFile = [BEHAV_IN taskBlock_resultsFileFndr(1).name];
        taskBlock_results = load(taskBlock_resultsFile);
        stimCondition = taskBlock_results.trialInfo(1).condition;
    end
    %stimIdx = find(strcmpi(TACS_CONDITIONS,stimCondition),1);

    % Figure out rest order
    % Commonly used reference string
    subDay = [subject '_' day];
    
    % preprocessing file
    preprocFile = [PROC_EEG_IN preprocPrefix subDay '.mat'];
    
    % Timing information
    if localComputer
        timingFile = [PROC_EEG_IN 'epochRejection.mat'];
    else
        timingFile = [PROC_EEG_IN subDay '_epochRejection.mat'];
    end
    
    %% CANONICAL PAC ANALYSIS
    % Loop seeds
    for seedIdx = 1:numSeedROIs_canonical
        seedRoi = canonical_SEED_ROIS{seedIdx};

        % Seed roi channels
        roiIdx = find(strcmp(seedRoi,ROI_DEF),1);
        roi_chanIdxs = ROI_DEF{roiIdx+1};

        % Loop pac pairs
        for pacPairIdx = 1:numPacPairs
            pacPair = PAC_PAIRS{pacPairIdx};


            canonicalPac_subjectFile = sprintf('%sseedPac_%s_%s_REST_%s_%s_%iiter_%s.mat',...
                PAC_OUTPUT,subject,stimCondition,seedRoi,pacPair,numIterations,analysisVersion);

            if exist(canonicalPac_subjectFile,'file')~=2

                clear canonicalPacStruct;
                % resting time series
                restDataStruct = load(preprocFile);
                restData_all = restDataStruct.restData;
                restName_all = restDataStruct.restNames;

                restTimingStruct = load(timingFile);
                restTiming_all = restTimingStruct.epochTiming;

                % loop through resting periods
                for restIdx = 1:numRestPeriods

                    % gather relevant information for this rest period
                    restName = restName_all{restIdx};
                    cond_data = restData_all{restIdx};
                    if isempty(cond_data)
                        continue;
                    end
                    restTiming_time = restTiming_all{1,restIdx};

                    taskCondIdxs = find(strcmpi(restName_all,restName));
                    firstOrSecondIter = find(restIdx == taskCondIdxs,1);

                    % Name the task condition more accurately
                    if length(find(taskCondIdxs)) == 1
                        taskCondition = restName;
                    elseif firstOrSecondIter == 1
                        taskCondition = ['first_' restName];
                    elseif firstOrSecondIter == 2
                        taskCondition = ['second_' restName];
                    else
                        error('should not reach here');
                    end

                    % Results from PAC analysis to be saved
                    % Raw values
                    PAC_allCond_allChan = NaN(numDataChannels,numMethods);
                    ampByPhase_allCond_allChan = NaN(numDataChannels,numPhaseBins);
                    % After iteration
                    PACz_allCond_allChan = NaN(numDataChannels,numMethods);
                    iterPAC_allCond_allChan = NaN(numDataChannels,numIterations,numMethods);

                    fprintf('\tfor %s %s %s %s %s\n',subject,stimCondition,taskCondition,seedRoi,pacPair);

                    % Seed ROI - average across channels (if applicable)
                    seed_eegData = squeeze(mean(cond_data(roi_chanIdxs,:),1));

                    % PAC pair information
                    pacPairParts = regexp(pacPair,'_','split');
                    lowFreqName  = pacPairParts{1};
                    highFreqName = pacPairParts{2};

                    % Frequencies range
                    lowFreqRange  = eval(lowFreqName);
                    highFreqRange = eval(highFreqName);

                    % calculate phase of low frequency
                    % Hilbert transform by default does a
                    % mirror flip of the data
                    [lowFreq_phase,~] = hilbertTransform_JR(...
                        sampleRate,lowFreqRange,seed_eegData);

                    % Loop data channels
                    for dataChanIdx = 1:numDataChannels
                        chanIdx = dataChannels(dataChanIdx);

                        fprintf('\t\t%i of %i channels\n',dataChanIdx,numDataChannels);
                        
                        target_eegData = squeeze(cond_data(chanIdx,:,:));

                        % Calculate high frequency amplitude
                        [~,highFreq_power] = hilbertTransform_JR(...
                            sampleRate,highFreqRange,target_eegData);

                        % Calculate low freuqency phase of high
                        % frequency amplitude
                        [highFreq_lowPhase,~] = hilbertTransform_JR(...
                            sampleRate,lowFreqRange,highFreq_power);

                        % Run PAC calculation - for all three methods
                        [PACz,PAC,iteratedPAC,ampByPhase] = PACz_all3_JR(...
                            lowFreq_phase,highFreq_power,highFreq_lowPhase,pacOptions);

                        % Store output of iterated PAC
                        PAC_allCond_allChan(dataChanIdx,:) = PAC;
                        PACz_allCond_allChan(dataChanIdx,:) = PACz;
                        ampByPhase_allCond_allChan(dataChanIdx,:) = ampByPhase;
                        iterPAC_allCond_allChan(dataChanIdx,:,:) = iteratedPAC;

                    end % loop channels

                    % Store data from this rest period
                    canonicalPacStruct.(taskCondition) = struct(...
                        'PAC',PAC_allCond_allChan,...
                        'PACz',PACz_allCond_allChan,...
                        'ampByPhase',ampByPhase_allCond_allChan,...
                        'iteratedPAC',iterPAC_allCond_allChan);

                end % loop rest periods

                % Save file out
                save(canonicalPac_subjectFile,'-struct','canonicalPacStruct');
            end % if file exists
        end % loop pac pair
    end % loop seeds


    %% EXHAUSTIVE PAC ANALYSIS
    % Loop seeds
    for pairRoiIdx = 1:numPairsROI_exhaustive
        seed_target = exhaustive_SEED_TARGET{pairRoiIdx};
        seed_target_parts = regexp(seed_target,'_','split');
        seedRoi   = seed_target_parts{1};
        targetRoi = seed_target_parts{2};


        % Does this file exist
        exhaustivePac_subjectFile = sprintf('%sexhPac_%s_%s_REST_seed%s_target%s_%iiter_%s.mat',...
            PAC_OUTPUT,subject,stimCondition,seedRoi,targetRoi,numIterations,analysisVersion);

        if exist(exhaustivePac_subjectFile,'file')~=2

            fprintf('Exhaustive PAC for %s %s %s\n',subject,seedRoi,targetRoi);

            % Seed roi channels
            seed_roiIdx = find(strcmp(seedRoi,ROI_NAMES),1);
            seedRoi_chanIdxs = ROI_DEF{seed_roiIdx*2};
            % Target roi channels
            target_roiIdx = find(strcmp(targetRoi,ROI_NAMES),1);
            targetRoi_chanIdxs = ROI_DEF{target_roiIdx*2};

            allCond_Files = cell(1,numRestPeriods);
            allTask_conditions = cell(1,numRestPeriods);

            % resting time series
            restDataStruct = load(preprocFile);
            restData_all = restDataStruct.restData;
            restName_all = restDataStruct.restNames;


            % loop through resting periods
            for restIdx = 1:numRestPeriods


                % gather relevant information for this rest period
                restName = restName_all{restIdx};

                taskCondIdxs = find(strcmpi(restName_all,restName));
                firstOrSecondIter = find(restIdx == taskCondIdxs,1);

                % Name the task condition more accurately
                if length(find(taskCondIdxs)) == 1
                    taskCondition = restName;
                elseif firstOrSecondIter == 1
                    taskCondition = ['first_' restName];
                elseif firstOrSecondIter == 2
                    taskCondition = ['second_' restName];
                else
                    error('should not reach here');
                end


                % Does this file exist
                cond_exhaustivePac_subjectFile = sprintf('%sexhPac_%s_%s_%s_REST_seed%s_target%s_%iiter_%s.mat',...
                    PAC_OUTPUT,subject,stimCondition,taskCondition,seedRoi,targetRoi,numIterations,analysisVersion);
                allCond_Files{restIdx} = cond_exhaustivePac_subjectFile;
                allTask_conditions{restIdx} = taskCondition;
                cond_data = restData_all{restIdx};

                if exist(cond_exhaustivePac_subjectFile,'file')~=2 && (~isempty(cond_data))

                    fprintf('Exhaustive PAC for %s %s %s\n',subject,seedRoi,targetRoi);

                    clear cond_exhaustivePacStruct;

                    % Results from PAC analysis to be saved
                    PAC_allCond_allChan = NaN(numLowFreq,numHighFreq,numMethods);
                    ampByPhase_allCond_allChan = NaN(numLowFreq,numHighFreq,numPhaseBins);
                    % After iteration
                    PACz_allCond_allChan = NaN(numLowFreq,numHighFreq,numMethods);
                    iterPAC_allCond_allChan = NaN(numLowFreq,numHighFreq,numIterations,numMethods);

                    % Seed ROI - average signal
                    seed_eegData = squeeze(mean(cond_data(seedRoi_chanIdxs,:),1));
                    % Target ROI - average signal
                    target_eegData = squeeze(mean(cond_data(targetRoi_chanIdxs,:),1));
                    
                    count = 0;

                    % Loop low frequencies for the seed
                    for lowFreqIdx = 1:numLowFreq
                        lowFreq = LOW_FREQ(lowFreqIdx);

                        % calculate low frequency phase
                        [lowFreq_phase,~]=phaseAmp_calculator(sampleRate,lowFreq,seed_eegData,numCycles_low);

                        % Loop high freq for the target
                        for highFreqIdx = 1:numHighFreq
                            highFreq = HIGH_FREQ(highFreqIdx);
                            
                            count = count + 1;
                            fprintf('\t\t%i of %i freq pairs\n',count,(numHighFreq*numLowFreq));

                            % Calculate high frequency amplitude
                            [~,highFreq_amp]=phaseAmp_calculator(...
                                sampleRate,highFreq,target_eegData,numCycles_high);

                            % Calculate low frequency phase of high frequency amplitude
                            [highFreq_lowPhase,~] = phaseAmp_calculator(...
                                sampleRate,lowFreq,highFreq_amp,numCycles_low);

                            % Calculate PAC
                            [PACz,PAC,iteratedPAC,ampByPhase] = PACz_all3_JR(...
                                lowFreq_phase,highFreq_amp,highFreq_lowPhase,pacOptions);

                            % Store output of iterated PAC
                            PAC_allCond_allChan(lowFreqIdx,highFreqIdx,:) = PAC;
                            PACz_allCond_allChan(lowFreqIdx,highFreqIdx,:) = PACz;
                            ampByPhase_allCond_allChan(lowFreqIdx,highFreqIdx,:) = ampByPhase;
                            iterPAC_allCond_allChan(lowFreqIdx,highFreqIdx,:,:) = iteratedPAC;

                        end % loop high freq
                    end % loop low freq
                    % Store data out for this rest period

                    cond_exhaustivePacStruct.(taskCondition) = struct(...
                        'PAC',PAC_allCond_allChan,...
                        'PACz',PACz_allCond_allChan,...
                        'ampByPhase',ampByPhase_allCond_allChan,...
                        'iteratedPAC',iterPAC_allCond_allChan);

                    save(cond_exhaustivePac_subjectFile,'-struct','cond_exhaustivePacStruct');

                end % loop conditions

            end % if file exists



            % resting data names
            clear exhaustivePacStruct;
            for restIdx = 1:numRestPeriods
                taskCondition = allTask_conditions{restIdx};

                cond_exhaustivePac_subjectFile = allCond_Files{restIdx};
                if exist(cond_exhaustivePac_subjectFile,'file')==2
                    cond_exhStruct = load(cond_exhaustivePac_subjectFile);
                    exhaustivePacStruct.(taskCondition) = cond_exhStruct;
                end
            end
            exhaustivePacStruct.order = {allTask_conditions};
            % Save data out
            save(exhaustivePac_subjectFile,'-struct','exhaustivePacStruct');
        end % save file out
    end % loop ROI pair for exhaustive PAC
end % loop days
end % end of function