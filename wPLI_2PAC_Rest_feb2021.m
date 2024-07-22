
analysisVersion = 'feb2021';

% Subjects
SUBJECTS = {...
    '01','02','03','04',...
    '05','06','07','08',...
    '09','10','11','12',...
    '14','15','16',... '13' non-compliant subject
    '17','18','19','20',...
    '21','22','23','24'};
numSub = length(SUBJECTS);


% Version of data
preprocVersion = 'apr2019';

% root
ROOT = '/Users/justinri/Dropbox (Frohlich Lab)/';

% Toolboxes
CODE = [ROOT 'Frohlich Lab Team Folder/Codebase/CodeJustin/'];

% Only load this locally
EEGLAB_TOOLBOX = [CODE 'Toolboxes/eeglab14_1_2b/'];
addpath(EEGLAB_TOOLBOX);
eeglab('rebuild');

RIDDLER_TOOLBOX = [CODE 'Riddler_Toolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));

% Data paths
DATA = [ROOT 'HumanStudies/2PAC/'];
RAW_BEHAV = [DATA 'RawBehavior/'];
PROC_EEG = [DATA 'ProcEEG/'];
GROUP_EEG = [DATA 'GroupEEG/'];
mkdir_JR(GROUP_EEG);

GROUP_CONN = [GROUP_EEG 'wPLI_Rest_' analysisVersion '/'];
mkdir_JR(GROUP_CONN);
SINGLE_CONN = sprintf('%sIndividual_AnalysisEEG/wPLI_Rest_%s/',DATA,analysisVersion);
mkdir_JR(SINGLE_CONN);


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
    'aParOcc',[37 42 53 54 61 78 79 86 87 93],...
    'focLM1', 36}; % derived from theta-gamma coupling
numROIs = length(ROI_DEF) / 2;
ROI_NAMES = ROI_DEF(1:2:length(ROI_DEF));

%% Hilbert transform in canonical bands
delta  = [ 2  3];
theta  = [ 4  7];
beta   = [15 25];
lowBeta= [10 20];
gamma  = [35 58];
% Connectivty analysis - interested in left and right motor - beta freq
CONN_FREQS = {'beta','delta','theta','gamma'};
numConnFreqs = length(CONN_FREQS);
conn_SEED_ROIS = {'lM1','rM1','cPFC','aPFC','rParOcc','focLM1'};
numSeedROIs_conn = length(conn_SEED_ROIS);

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



%% Loop through participants and calculate connectivity
for subIdx = 18:numSub
    subject =['sub' SUBJECTS{subIdx}];

    % output file for single subject level
    CONN_OUTPUT = [SINGLE_CONN subject '/'];
    mkdir_JR(CONN_OUTPUT);

    SUB_PROC_EEG  = [PROC_EEG subject '/'];

    SUB_RAW_BEHAV = [RAW_BEHAV subject '/'];

    % Loop through days
    for dayIdx = 1:numDays
        day = STIM_DAYS{dayIdx};

        % Subject and day raw behavior
        BEHAV_IN = [SUB_RAW_BEHAV day '/'];

        % Preprocessing directory
        PROC_EEG_IN = [SUB_PROC_EEG day '_Rest_' preprocVersion '/'];

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
        timingFile = [PROC_EEG_IN 'epochRejection.mat'];

        %% Seed-based connectivity analysis
        % Loop seeds
        for seedIdx = 1:numSeedROIs_conn
            seedRoi = conn_SEED_ROIS{seedIdx};

            % Seed roi channels
            roiIdx = find(strcmp(seedRoi,ROI_DEF),1);
            roi_chanIdxs = ROI_DEF{roiIdx+1};

            % Loop CONN pairs
            for connFreqIdx = 1:numConnFreqs
                connFreq = CONN_FREQS{connFreqIdx};


                seedConn_subjectFile = sprintf('%sseedConn_%s_%s_REST_%s_%s_%s.mat',...
                    CONN_OUTPUT,subject,stimCondition,seedRoi,connFreq,analysisVersion);

                if exist(seedConn_subjectFile,'file')~=2

                    clear seedConnStruct;
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

                        % Results from CONN analysis to be saved
                        % Raw values
                        CONN_allCond_allChan = NaN(numDataChannels,1);
                        fprintf('\tfor %s %s %s %s %s\n',subject,stimCondition,taskCondition,seedRoi,connFreq);

                        % Seed ROI - average across channels (if applicable)
                        if length(roi_chanIdxs) > 1
                            seed_eegData = squeeze(mean(cond_data(roi_chanIdxs,:),1));
                        else
                            seed_eegData = squeeze(cond_data(roi_chanIdxs,:));
                        end

                        % CONN pair information
                        % Frequencies range
                        freqRange  = eval(connFreq);

                        % calculate phase of low frequency
                        % Hilbert transform by default does a
                        % mirror flip of the data
                        [seedConn_phase,seedConn_amp] = hilbertTransform_JR(...
                            sampleRate,freqRange,seed_eegData);
                        seedConn_analytic = ampPhase_to_analyticSignal(seedConn_amp,seedConn_phase);

                        % Loop data channels
                        for dataChanIdx = 1:numDataChannels
                            chanIdx = dataChannels(dataChanIdx);

                            %fprintf('\t\t%i of %i channels\n',dataChanIdx,numDataChannels);
                            target_eegData = squeeze(cond_data(chanIdx,:));

                            % Calculate high frequency amplitude
                            [targetConn_phase,targetConn_amp] = hilbertTransform_JR(...
                                sampleRate,freqRange,target_eegData);
                            targetConn_analytic = ampPhase_to_analyticSignal(targetConn_amp,targetConn_phase);

                            % Run CONN calculation - for all three methods
                            [~,wPLI,~] = phaseLagIndex_JR(...
                                seedConn_analytic,targetConn_analytic);

                            % Store output of iterated PAC
                            CONN_allCond_allChan(dataChanIdx) = wPLI;

                        end % loop channels

                        % Store data from this rest period
                        seedConnStruct.(taskCondition) = struct(...
                            'wPLI',CONN_allCond_allChan);

                    end % loop rest periods

                    % Save file out
                    save(seedConn_subjectFile,'-struct','seedConnStruct');
                end % if file exists
            end % loop pac pair
        end % loop seeds
    end % loop days
end % loop subjects