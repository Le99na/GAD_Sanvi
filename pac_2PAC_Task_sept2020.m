function pac_2PAC_Task_sept2020(subIdx,varargin)

% Important options
highDefinitionFLAG = 0;
diffAnalysisWindowFLAG = 1;

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
preprocVersion = 'oct2018';

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
    PROC_EEG = [DATA 'ProcEEG/'];
    GROUP_EEG = [DATA 'GroupEEG/'];
    mkdir_JR(GROUP_EEG);
    GROUP_PAC = [GROUP_EEG 'PAC_Task_' analysisVersion '/'];
    mkdir_JR(GROUP_PAC);
    SINGLE_PAC = sprintf('%sIndividual_AnalysisEEG/PAC_Task_%s/',DATA,analysisVersion);
    mkdir_JR(SINGLE_PAC);

    %% Loop through participants and calculate PAC
    subject =['sub' SUBJECTS{subIdx}];
    BASELINE_EEG = [PROC_EEG subject '/Baseline_' preprocVersion  '/'];

    % output file for single subject level
    PAC_OUTPUT = [SINGLE_PAC subject '/'];
    mkdir_JR(PAC_OUTPUT);
    
else
    
    %% LONGLEAF CLUSTER
    % NAS Root
    NAS_ROOT = '/nas/longleaf/home/justinri/';
    % Pine Root
    PINE_ROOT = '/pine/scr/j/u/justinri/2PAC/';

    % Code base on cluster
    CODE = [NAS_ROOT 'CodeJustin/'];
    
    % Data paths
    BASELINE_EEG = [PINE_ROOT 'ProcEEG_Baseline/'];

    % output file for single subject level
    PAC_OUTPUT =[PINE_ROOT 'PAC_Baseline/'];
     
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
chanlocFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_channelLocations.mat'];
chanlocsStruct = load(chanlocFile);
chanlocs = chanlocsStruct.chanlocs;
dataChannelsFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_dataChannels_idxs.mat'];
dataChannelsStruct = load(dataChannelsFile);
dataChannels = dataChannelsStruct.dataChannels;
numDataChannels = length(dataChannels);
data_chanlocs = chanlocs(dataChannels);

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
ROI_NAMES = ROI_DEF(1:2:numROIs);

%% Version 1 - Hilbert transform in canonical bands
delta  = [ 2  3];
theta  = [ 4  7];
beta   = [15 25];
lowBeta= [10 20];
gamma  = [35 58];
% Two specific hypotheses
PAC_PAIRS = {'delta_beta','theta_gamma'};
numPacPairs = length(PAC_PAIRS);
canonical_SEED_ROIS = {'cPFC','aPFC','rPFC'};
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

% Analysis window
if diffAnalysisWindowFLAG
    %analysisWindow = [0.2 1.4];
    lowFreqWindow = {...
        [0.2 1.4],...
        [0.6 1.2]};
    lowCut = 3.999;
    lowCutStr = 'lowCut_';
else
    analysisWindow = [0.4 1.6];
    lowCutStr = '';
end

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

        canonicalPac_subjectFile = sprintf('%s%sseedPac_%s_%s_%s_%iiter_%s.mat',...
            PAC_OUTPUT,lowCutStr,subject,seedRoi,pacPair,numIterations,analysisVersion);

        if exist(canonicalPac_subjectFile,'file')~=2

            % Results from PAC analysis to be saved
            % Raw values
            PAC_allCond_allChan = NaN(numCond,numDataChannels,numMethods);
            ampByPhase_allCond_allChan = NaN(numCond,numDataChannels,numPhaseBins);
            % After iteration
            PACz_allCond_allChan = NaN(numCond,numDataChannels,numMethods);
            iterPAC_allCond_allChan = NaN(numCond,numDataChannels,numIterations,numMethods);

            % Loop through each in condition
            for condIdx = 1:numCond
                cond = COND{condIdx};

                % Load preprocessed epoched data
                preprocFile = [BASELINE_EEG 'final_ica_' cond '_labeled_preproc_2PAC_' subject '_Baseline'];

                % Load condition data
                if localComputer
                    EEG = pop_loadset([preprocFile '.set']);
                else
                    EEG = load([preprocFile '.mat']);
                end
                cond_data = EEG.data;
                numEpochs = size(cond_data,3);

                fprintf('\tCanonical PAC: %s %s %s %s\n',subject,cond,seedRoi,pacPair);

                % Seed ROI - average across channels (if applicable)
                seed_eegData = squeeze(mean(cond_data(roi_chanIdxs,:,:),1));

                % PAC pair information
                pacPairParts = regexp(pacPair,'_','split');
                lowFreqName  = pacPairParts{1};
                highFreqName = pacPairParts{2};

                % Frequencies range
                lowFreqRange  = eval(lowFreqName);
                highFreqRange = eval(highFreqName);
                
                if diffAnalysisWindowFLAG
                    if lowFreqRange(1) > lowCut
                        analysisWindow = lowFreqWindow{2};
                    else
                        analysisWindow = lowFreqWindow{1};
                    end
                end
                
                window_startIdx = find((realTime+0.0001) >= analysisWindow(1),1);
                window_stopIdx  = find((realTime+0.0001) >= analysisWindow(2),1)-1;
                windowIdxs = window_startIdx:window_stopIdx;
                numWindowPts = length(windowIdxs);
                
                % calculate phase of low frequency
                % Hilbert transform by default does a
                % mirror flip of the data
                all_lowFreq_phase = NaN(numTimePoints,numEpochs);
                for epochIdx = 1:numEpochs
                    [lowFreq_phase,~] = hilbertTransform_JR(...
                        sampleRate,lowFreqRange,seed_eegData(:,epochIdx));
                    all_lowFreq_phase(:,epochIdx) = lowFreq_phase;
                end
                win_lowFreq_phase = all_lowFreq_phase(windowIdxs,:);
                lowFreq_windowPhase = reshape(win_lowFreq_phase,1,numWindowPts * numEpochs);

                % Loop data channels
                for dataChanIdx = 1:numDataChannels
                    chanIdx = dataChannels(dataChanIdx);
                    fprintf('\tchan %i of %i\n',dataChanIdx,numDataChannels);
                    
                    target_eegData = squeeze(cond_data(chanIdx,:,:));

                    all_highFreq_power = NaN(numTimePoints,numEpochs);
                    all_highFreq_lowPhase = NaN(numTimePoints,numEpochs);
                    for epochIdx = 1:numEpochs

                        % Calculate high frequency amplitude
                        [~,highFreq_power] = hilbertTransform_JR(...
                            sampleRate,highFreqRange,target_eegData(:,epochIdx));
                        all_highFreq_power(:,epochIdx) = highFreq_power;

                        % Calculate low freuqency phase of high
                        % frequency amplitude
                        [highFreq_lowPhase,~] = hilbertTransform_JR(...
                            sampleRate,lowFreqRange,highFreq_power);
                        all_highFreq_lowPhase(:,epochIdx) = highFreq_lowPhase;

                    end
                    win_highFreq_power = all_highFreq_power(windowIdxs,:);
                    highFreq_windowPower = reshape(win_highFreq_power,1,numWindowPts * numEpochs);

                    win_highFreq_lowPhase = all_highFreq_lowPhase(windowIdxs,:);
                    highFreq_windowLowPhase = reshape(win_highFreq_lowPhase,1,numWindowPts * numEpochs);

                    % Run PAC calculation - for all three methods
                    [PACz,PAC,iteratedPAC,ampByPhase] = PACz_all3_JR(...
                        lowFreq_windowPhase,highFreq_windowPower,highFreq_windowLowPhase,pacOptions);

                    % Store output of iterated PAC
                    PAC_allCond_allChan(condIdx,dataChanIdx,:) = PAC;
                    PACz_allCond_allChan(condIdx,dataChanIdx,:) = PACz;
                    ampByPhase_allCond_allChan(condIdx,dataChanIdx,:) = ampByPhase;
                    iterPAC_allCond_allChan(condIdx,dataChanIdx,:,:) = iteratedPAC;

                end % loop channels
            end % loop conditions
            % Save data out
            canonicalPacStruct = struct(...
                'PAC',PAC_allCond_allChan,...
                'PACz',PACz_allCond_allChan,...
                'ampByPhase',ampByPhase_allCond_allChan,...
                'iteratedPAC',iterPAC_allCond_allChan);
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

    % Seed roi channels
    seed_roiIdx = find(strcmp(seedRoi,ROI_NAMES),1);
    seedRoi_chanIdxs = ROI_DEF{seed_roiIdx*2};
    % Target roi channels
    target_roiIdx = find(strcmp(targetRoi,ROI_NAMES),1);
    targetRoi_chanIdxs = ROI_DEF{target_roiIdx*2};

    % Does this file exist
    exhaustivePac_subjectFile = sprintf('%s%sexhPac_%s_seed%s_target%s_%iiter_%inumFreq_%s.mat',...
        PAC_OUTPUT,lowCutStr,subject,seedRoi,targetRoi,numIterations,numFreqs,analysisVersion);

    if exist(exhaustivePac_subjectFile,'file')~=2
        
        fprintf('Exhaustive PAC for %s %s %s\n',subject,seedRoi,targetRoi);

        % Results from PAC analysis to be saved
        PAC_allCond_allChan = NaN(numCond,numLowFreq,numHighFreq,numMethods);
        ampByPhase_allCond_allChan = NaN(numCond,numLowFreq,numHighFreq,numPhaseBins);
        % After iteration
        PACz_allCond_allChan = NaN(numCond,numLowFreq,numHighFreq,numMethods);
        iterPAC_allCond_allChan = NaN(numCond,numLowFreq,numHighFreq,numIterations,numMethods);

        % Loop through each in condition
        for condIdx = 1:numCond
            cond = COND{condIdx};

            % Load preprocessed epoched data
            preprocFile = [BASELINE_EEG 'final_ica_' cond '_labeled_preproc_2PAC_' subject '_Baseline'];

            % Load condition data
            if localComputer
                EEG = pop_loadset([preprocFile '.set']);
            else
                EEG = load([preprocFile '.mat']);
            end
            cond_data = EEG.data;
            numEpochs = size(cond_data,3);

            % Seed ROI - average signal
            seed_eegData = squeeze(mean(cond_data(seedRoi_chanIdxs,:,:),1));
            % Target ROI - average signal
            target_eegData = squeeze(mean(cond_data(targetRoi_chanIdxs,:,:),1));

            count = 0;
            
            % Loop low frequencies for the seed
            for lowFreqIdx = 1:numLowFreq
                lowFreq = LOW_FREQ(lowFreqIdx);

                
                if diffAnalysisWindowFLAG
                    if lowFreqRange(1) > lowCut
                        analysisWindow = lowFreqWindow{2};
                    else
                        analysisWindow = lowFreqWindow{1};
                    end
                end
                
                window_startIdx = find((realTime+0.0001) >= analysisWindow(1),1);
                window_stopIdx  = find((realTime+0.0001) >= analysisWindow(2),1)-1;
                windowIdxs = window_startIdx:window_stopIdx;
                numWindowPts = length(windowIdxs);
                
                % calculate low frequency phase
                [lowFreq_phase,~]=phaseAmp_mirrorCalc(sampleRate,lowFreq,seed_eegData,numCycles_low);
                win_lowFreq_phase = lowFreq_phase(windowIdxs,:);
                lowFreq_windowPhase = reshape(win_lowFreq_phase,1,numWindowPts * numEpochs);

                % Loop high freq for the target
                for highFreqIdx = 1:numHighFreq
                    highFreq = HIGH_FREQ(highFreqIdx);

                    % Calculate high frequency amplitude
                    [~,highFreq_amp]=phaseAmp_mirrorCalc(...
                        sampleRate,highFreq,target_eegData,numCycles_high);
                    win_highFreq_power = highFreq_amp(windowIdxs,:);
                    highFreq_windowPower = reshape(win_highFreq_power,1,numWindowPts * numEpochs);

                    % Calculate low frequency phase of high frequency amplitude
                    [highFreq_lowPhase,~] = phaseAmp_mirrorCalc(...
                        sampleRate,lowFreq,highFreq_amp,numCycles_low);
                    win_highFreq_lowPhase = highFreq_lowPhase(windowIdxs,:);
                    highFreq_windowLowPhase = reshape(win_highFreq_lowPhase,1,numWindowPts * numEpochs);

                    count = count + 1;
                    fprintf('\t%i of %i freq pairs\n',count,numLowFreq * numHighFreq);
                    
                    % Calculate PAC
                    [PACz,PAC,iteratedPAC,ampByPhase] = PACz_all3_JR(...
                        lowFreq_windowPhase,highFreq_windowPower,highFreq_windowLowPhase,pacOptions);

                    % Store output of iterated PAC
                    PAC_allCond_allChan(condIdx,lowFreqIdx,highFreqIdx,:) = PAC;
                    PACz_allCond_allChan(condIdx,lowFreqIdx,highFreqIdx,:) = PACz;
                    ampByPhase_allCond_allChan(condIdx,lowFreqIdx,highFreqIdx,:) = ampByPhase;
                    iterPAC_allCond_allChan(condIdx,lowFreqIdx,highFreqIdx,:,:) = iteratedPAC;

                end % loop high freq
            end % loop low freq
        end % loop conditions
        % Save data out
        exhaustivePacStruct = struct(...
            'PAC',PAC_allCond_allChan,...
            'PACz',PACz_allCond_allChan,...
            'ampByPhase',ampByPhase_allCond_allChan,...
            'iteratedPAC',iterPAC_allCond_allChan);
        save(exhaustivePac_subjectFile,'-struct','exhaustivePacStruct');
    end % if file exists
end % loop ROI pair for exhaustive PAC

end % function ends