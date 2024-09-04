% wipe out all previous data
clear all; clc; close all;

% Preprocessing version
preprocVersion = '073024_3';

% Use ICLabel or not
useIClabel = 1;
skipManualIC = 0;

%% Step 0: get organized
% Paths to data
root = 'C:\Users\sanvi\Documents\GAD EEG ERP Analysis';
addpath(genpath(root))
ClassObj = NNC_path();

%generate cell array of available participants
grps_strct = dir(ClassObj.DataRoot_GAD_anxious_original);
grp_nams = {grps_strct.name};

longerThan3 = cellfun(@(x) length(x) > 3, grp_nams);
GROUPS = grp_nams(find(longerThan3));
grps_2 = {};
for i=1:length(GROUPS)
    grps_2{i} = [GROUPS{i} '\EEG\'];
end

numGroups = length(grps_2);

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

% Start EEGlab
eeglab

%% Information about the resting-state period
% Information about the filtering
downSampleRate  = 200;
highpassfreq    = 1;
lowpassfreq     = 40;

% for global average rereference give it an empty bracket
referenceElectrodes = [];

% display available data
grps_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompt user for ptp to process
groupIdx = input('Please enter a number for groupIdx: ');

% Name of the group which is used in folder organization
group = grps_2{groupIdx};
grp_name_lng = (split(group, '\'));
grp_nam = grp_name_lng{1};

% Create a folder system for the processed data
PROC_GROUP = [ClassObj.DataRoot_GAD_anxious_processed group '/'];
if exist(PROC_GROUP,'dir')~=7
    mkdir(PROC_GROUP);
end

% Get all task recordings from each subjects
lst_grp_tot = ClassObj.getRecs([ClassObj.DataRoot_GAD_anxious_original group]);

%look into task eeg file
foundFilename = lst_grp_tot{1};
rec = grp_nam;
path_fl_raw = [ClassObj.DataRoot_GAD_anxious_original group foundFilename];

PROC_GROUP = [PROC_GROUP rec '/'];
if exist(PROC_GROUP,'dir')~=7
    mkdir(PROC_GROUP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: Import the data

% import the raw EEG file
EEG = pop_mffimport(path_fl_raw,{'code'});
photodiodeTimeseries = EEG.data(130, :);

% define file for data saving
EEG.filename = [rec '_preVer_' preprocVersion '.set'];
eegFile = [PROC_GROUP EEG.filename(1:end-4)];

% Butterworth high pass filter 1Hz
EEG.data = butterFilt(EEG.data',EEG.srate,highpassfreq,'high')';

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4) '_filtHigh.set'],'filepath',PROC_GROUP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save photodiode data
peripheralsTime = EEG.times;
photodiodeFile = [EEG.filepath '\photo_GAD_task.mat'];
rawEventStruct = struct(...
    'photodiode',photodiodeTimeseries,...
    'time',peripheralsTime,...
    'events',EEG.event);
save(photodiodeFile,'-struct','rawEventStruct');

% Delete trigger channel
EEG = pop_select( EEG, 'nochannel',{'Trigger','PhotoDiode'});

% Saves times we cut
EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4) '_delPhoto.set'],'filepath',PROC_GROUP);

ALLEEG{1} = EEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity check - create spectrogram of raw data
% Setting paramter for psd sanity checks
Flim = [1 80];                          % (Hz)
taper = 1;                              % Parameter for Tukey Window 1=cosine
foi         = Flim(1):0.1:Flim(2);      % (Hz) frequencies of interest
Fs          = EEG.srate;
winSize     = 2*Fs;                     % (sample) window size for computing spectra
win         = tukeywin(winSize, taper); % the Tukey window
lag         = round(Fs/8);              % (sample) lag between windows
overlap     = winSize - lag;            % (sample) overlap between windows

% -- computing Welch's Power Spectral Density Estimate
disp("Computing PSD. Perhaps slow...")
tic
[pxx,a] = pwelch(detrend(EEG.data'),win, overlap, foi, Fs);
toc

% Plotting
figure
freqRange = foi >= 8 & foi <= 12;  % Frequency indices within 8-12 Hz range
maxPower = max(max(pxx(freqRange, :)));  % Maximum power in the range
cl=colorlines(plot(foi,pxx),eegcmap([],EEG.chanlocs));
labs = {EEG.chanlocs.labels};
[cl.Tag] = labs{:};
[cl.ButtonDownFcn] = deal(@(x,~) disp(x.Tag));
xlim(Flim)
ylim([0 maxPower+10])
ylabel('PSD [\mu V^2/Hz]')
xlabel('Frequency [Hz]')
title('Spectrum after cutting & high pass')
plotElectHead(EEG.chanlocs)

% Save png
filename_step1 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step1' , '.png'];
saveas(gcf, filename_step1, 'png');

filename_step1_1 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step1' , '.fig'];
savefig(gcf,filename_step1_1)

confirmation = input('Do you want to close all figures and continue? (y/n): ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('All figures are closed.');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 2: Filter and downsample

% Butterworth bandpass filter
EEG.data = butterFilt(EEG.data',EEG.srate,lowpassfreq,'low')';

% Run down-sampling
EEG = pop_resample(EEG,downSampleRate);

% Save the data
EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4) '_filtLow.set'],'filepath',EEG.filepath);

% Sanity check - create spectrogram of raw data
% Setting paramter for psd sanity checks
Flim = [1 80];                          % (Hz)
taper = 1;                              % Parameter for Tukey Window 1=cosine
foi         = Flim(1):0.1:Flim(2);      % (Hz) frequencies of interest
Fs          = EEG.srate;
winSize     = 2*Fs;                     % (sample) window size for computing spectra
win         = tukeywin(winSize, taper); % the Tukey window
lag         = round(Fs/8);              % (sample) lag between windows
overlap     = winSize - lag;            % (sample) overlap between windows

% -- computing Welch's Power Spectral Density Estimate
disp("Computing PSD. Perhaps slow...")
tic
[pxx,a] = pwelch(detrend(EEG.data'),win, overlap, foi, Fs);
toc

% Plotting
figure
freqRange = foi >= 8 & foi <= 12;  % Frequency indices within 8-12 Hz range
maxPower = max(max(pxx(freqRange, :)));  % Maximum power in the range
cl=colorlines(plot(foi,pxx),eegcmap([],EEG.chanlocs));
labs = {EEG.chanlocs.labels};
[cl.Tag] = labs{:};
[cl.ButtonDownFcn] = deal(@(x,~) disp(x.Tag));
xlim(Flim)
ylim([0 maxPower+10])
ylabel('PSD [\mu V^2/Hz]')
xlabel('Frequency [Hz]')
title('Spectrum after cutting & high pass')
plotElectHead(EEG.chanlocs)

% Save png
filename_step2 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step2' , '.png'];
saveas(gcf, filename_step2, 'png');

filename_step2_1 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step2' , '.fig'];
savefig(gcf,filename_step2_1)

confirmation = input('Do you want to close all figures and continue? (y/n): ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('All figures are closed.');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 3: Remove bad channels automized
%%%% Remove bad channels (ideally not more than 20%, i.e. 25 channels with 128 electrodes net)
% Adapt the events so that they are displayed in the scroll

% remove bad electrodes visually
EEG = csc_eeg_plotter(EEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prompt the user to check if any channels were rejected
channels_rejected = input('Did you reject any channels? (y/n): ', 's');

% If channels were rejected, run the interpolation and saving section
if strcmpi(channels_rejected, 'y')

    % Adapt the structure of the EEG
    bad_channels_idx = sort(EEG.hidden_channels);
    ch_labels = {EEG.chanlocs.labels};
    bad_ch_labels = ch_labels(bad_channels_idx);

    % update EEG structure
    EEG.bad_channels_labels = bad_ch_labels; % EEG.eBridge.Count
    EEG.bad_channels_index = sort([cellfun(@(x) str2double(x(2:end)), EEG.bad_channels_labels)]);
    EEG = epi_log(@pop_select, EEG, 'nochannel', bad_channels_idx);% remove bad channels (based on indices)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % visualize
    figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % check remaining electrodes (if many electrodes rejected in one area, might be problematic for interpolation and for results interpretation)
    figure; topoplot([],ALLEEG{1}.chanlocs(EEG.bad_channels_index), 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % which channels removed

    % backup
    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4) '_noBadChan1.set'],'filepath',EEG.filepath);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Interpolate rejected electrodes
    EEG = pop_interp(EEG,ALLEEG{1}.chanlocs, 'spherical');

    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ipol1.set'], 'filepath', EEG.filepath);

    % Sanity check - create spectrogram of raw data
    % Setting paramter for psd sanity checks
    Flim = [1 80];                          % (Hz)
    taper = 1;                              % Parameter for Tukey Window 1=cosine
    foi         = Flim(1):0.1:Flim(2);      % (Hz) frequencies of interest
    Fs          = EEG.srate;
    winSize     = 2*Fs;                     % (sample) window size for computing spectra
    win         = tukeywin(winSize, taper); % the Tukey window
    lag         = round(Fs/8);              % (sample) lag between windows
    overlap     = winSize - lag;            % (sample) overlap between windows

    % -- computing Welch's Power Spectral Density Estimate
    disp("Computing PSD. Ugh so slow...")
    tic
    [pxx,a] = pwelch(detrend(EEG.data'),win, overlap, foi, Fs);
    toc

    % Plotting
    figure
    freqRange = foi >= 8 & foi <= 12;  % Frequency indices within 8-12 Hz range
    maxPower = max(max(pxx(freqRange, :)));  % Maximum power in the range
    cl=colorlines(plot(foi,pxx),eegcmap([],EEG.chanlocs));
    labs = {EEG.chanlocs.labels};
    [cl.Tag] = labs{:};
    [cl.ButtonDownFcn] = deal(@(x,~) disp(x.Tag));
    xlim(Flim)
    ylim([0 maxPower+10])
    ylabel('PSD [\mu V^2/Hz]')
    xlabel('Frequency [Hz]')
    title('Spectrum after cutting & high pass')
    plotElectHead(EEG.chanlocs)

    % Save png
    filename_step3 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step3' , '.png'];
    saveas(gcf, filename_step3, 'png');

    filename_step3_1 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step3' , '.fig'];
    savefig(gcf,filename_step3_1)

    confirmation = input('Do you want to close all figures and continue? (y/n): ', 's');

    % Check the user's response
    if strcmpi(confirmation, 'y')
        close all;
        disp('All figures are closed.');
    else
        error('Processing discontinued.');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Second possibility for channel rejection (save only if you DID reject further channels)
    % Ask the user for input

    doSecondRejection = input('Do you want to perform the second channel rejection? (y/n): ', 's');

    if strcmpi(doSecondRejection, 'y')

        % remove bad electrodes visually
        EEG = csc_eeg_plotter(EEG);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Adapt the structure of the EEG
        bad_channels_idx = sort(EEG.hidden_channels);
        ch_labels={EEG.chanlocs.labels};
        bad_ch_labels=ch_labels(bad_channels_idx);

        % update EEG structure
        EEG.bad_channels_labels = [EEG.bad_channels_labels bad_ch_labels]; % EEG.eBridge.Count
        EEG.bad_channels_index = sort([cellfun(@(x) str2double(x(2:end)), EEG.bad_channels_labels)]);
        EEG = epi_log(@pop_select, EEG, 'nochannel', bad_channels_idx); %removing channels

        EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_noBadChan2.set'], 'filepath', EEG.filepath);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Second channel rejection (run only if channels rejected the second time)

        % Ask the user if they want to perform the second channel rejection

        % Interpolate rejected electrodes
        EEG = pop_interp(EEG,ALLEEG{1}.chanlocs, 'spherical');
        EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ipol2.set'], 'filepath', EEG.filepath);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Sanity check - create spectrogram of raw data
        % Setting paramter for psd sanity checks
        Flim = [1 40];                          % (Hz)
        taper = 1;                              % Parameter for Tukey Window 1=cosine
        foi         = Flim(1):0.1:Flim(2);      % (Hz) frequencies of interest
        Fs          = EEG.srate;
        winSize     = 2*Fs;                     % (sample) window size for computing spectra
        win         = tukeywin(winSize, taper); % the Tukey window
        lag         = round(Fs/8);              % (sample) lag between windows
        overlap     = winSize - lag;            % (sample) overlap between windows

        % -- computing Welch's Power Spectral Density Estimate
        disp("Computing PSD. Maybe sloooooooooooooow...")
        tic
        [pxx,~] = pwelch(detrend(EEG.data'),win, overlap, foi, Fs);
        toc

        % Plotting
        figure
        cl=colorlines(plot(foi,pxx),eegcmap([],EEG.chanlocs));
        labs = {EEG.chanlocs.labels};
        [cl.Tag] = labs{:};
        [cl.ButtonDownFcn] = deal(@(x,~) disp(x.Tag));
        xlim(Flim)
        ylim([0 60])
        ylabel('PSD [\mu V^2/Hz]')
        xlabel('Frequency [Hz]')
        title('After cutting+filter+interpolation')
        plotElectHead(EEG.chanlocs)

        % Save png and.fig for Step 2
        filename_step3_2 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step3.5' , '.png'];
        saveas(gcf, filename_step3_2, 'png');

        filename_step3_3 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step3.5' , '.fig'];
        savefig(gcf,filename_step3_3)


        confirmation = input('Do you want to close all figures and continue? (y/n): ', 's');

        % Check the user's response
        if strcmpi(confirmation, 'y')
            close all;
            disp('All figures are closed.');
        else
            error('Processing discontinued.');
        end

    else

        disp('Skipping second channel rejection.');

    end
else

    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4) '_noChanRej.set'],'filepath',EEG.filepath);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 4: Global average rereferencing

% Average re-referencing
referenceElectrodes = [];
EEG = pop_reref(EEG,referenceElectrodes);

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_GloAv.set'], 'filepath', EEG.filepath);

% Sanity check - create spectrogram of raw data
% Setting paramter for psd sanity checks
Flim = [1 40];                          % (Hz)
taper = 1;                              % Parameter for Tukey Window 1=cosine
foi         = Flim(1):0.1:Flim(2);      % (Hz) frequencies of interest
Fs          = EEG.srate;
winSize     = 2*Fs;                     % (sample) window size for computing spectra
win         = tukeywin(winSize, taper); % the Tukey window
lag         = round(Fs/8);              % (sample) lag between windows
overlap     = winSize - lag;            % (sample) overlap between windows

% -- computing Welch's Power Spectral Density Estimate
disp("Computing PSD. OMG so slow...")
tic
[pxx,~] = pwelch(detrend(EEG.data'),win, overlap, foi, Fs);
toc

% Plotting
figure
freqRange = foi >= 8 & foi <= 12;  % Frequency indices within 8-12 Hz range
maxPower = max(max(pxx(freqRange, :)));  % Maximum power in the range
cl=colorlines(plot(foi,pxx),eegcmap([],EEG.chanlocs));
labs = {EEG.chanlocs.labels};
[cl.Tag] = labs{:};
[cl.ButtonDownFcn] = deal(@(x,~) disp(x.Tag));
xlim(Flim)
ylim([0 maxPower+10])
ylabel('PSD [\mu V^2/Hz]')
xlabel('Frequency [Hz]')
title(['After cutting+filter+interpolation+average - Subject ', num2str(rec)]);
plotElectHead(EEG.chanlocs)

% Save png + .fig
filename_step4 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step4' , '.png'];
saveas(gcf, filename_step4, 'png');

filename_step4_1 = [EEG.filepath '\' EEG.filename(1:end-4) '_Step4' , '.fig'];
savefig(gcf,filename_step4_1)

confirmation = input('Do you want to close all figures and continue? (y/n): ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('All figures are closed.');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% use photodiode to repaire triggers!
b = load([EEG.filepath '\photo_GAD_task.mat']);

plot(b.time,b.photodiode)

photo_h = input('Check photodiode triggers and enter value for minimum trigger height: ', 's');
photo_h_num = str2num(photo_h);

idx_mss_all = ClassObj.getPhotodiode(b,EEG,photo_h_num);

idx_mss_all
confirmation = input('Check validity of missing triggers. Continue? (y/n) ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('Continuing...');
else
    error('Processing discontinued.');
end

%idx_mss = flip(idx_mss_all)

event_new = EEG.event;

for i = 1:length(idx_mss_all)
    idx = idx_mss_all(i);
    event_new = [event_new(1:idx-1), struct('description',[],'begintime',[],'classid',[],'code',[],'duration',[],'label',[],'relativebegintime',[],'sourcedevice',[],'latency',[],'type',[],'mffkeys',[],'mffkey_gidx',[],'mffkey_cidx',[],'mffkeysbackup',[]), event_new(idx:end)];
end

openvar('event_new')
confirmation = input('Check empty lines are inserted at stimulusstart or stimulusend. Continue? (y/n) ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    disp('Continuing...');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%% check in variable event_new  if empty lines are inserted in correct location
%%%%%%%%%%%%%% (should only be in location of stimulusstart or stimulusend)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in trigger labels
path_behav = [ClassObj.DataRoot_GAD_anxious_original grp_nam '\' 'Behavioral\'];
lst_trggr = dir(path_behav);
trggr_nams = {lst_trggr.name};
trggr_fl = trggr_nams{find(contains(trggr_nams,'block_11_trigger_labels'))};
addpath(path_behav)
lbls = readtable(trggr_fl,'ReadVariableNames', false);
lst_xlsx=[lbls.Var1;lbls.Var2;lbls.Var3;lbls.Var4;lbls.Var5;lbls.Var6;lbls.Var7;lbls.Var8;lbls.Var9;lbls.Var10;lbls.Var11;lbls.Var12];
lst_xlsx_def=lst_xlsx(find(cellfun(@(x) length(x)>2, lst_xlsx)));

lst_trr_c = {event_new.latency};
%transformation of latencies into ms
lst_lat_ms = cellfun(@(x)(x/EEG.srate)*1000, lst_trr_c,'UniformOutput', false);

idx_opt= ClassObj.getBinaries(lst_lat_ms,lst_xlsx);

idx_opt
confirmation = input('Check validity of missing triggers. Continue? (y/n) ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('Continuing...');
else
    error('Processing discontinued.');
end

rev_idx_opt = flip(idx_opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(idx_opt)  %careful:sometimes you have to manually adjust EEG.event and delete columns
    if rev_idx_opt(j)==1
        % Shift elements below the insertion point down by one
        event_new = [struct('description',[],'begintime',[],'classid',[],'code',[],'duration',[],'label',[],'relativebegintime',[],'sourcedevice',[],'latency',[],'type',[],'mffkeys',[],'mffkey_gidx',[],'mffkey_cidx',[],'mffkeysbackup',[]), event_new(1:end)];
    else
        num = num2str(rev_idx_opt(j));
        event_new = [event_new(1:rev_idx_opt(j)-1), struct('description',[],'begintime',[],'classid',[],'code',[],'duration',[],'label',[],'relativebegintime',[],'sourcedevice',[],'latency',[],'type',[],'mffkeys',[],'mffkey_gidx',[],'mffkey_cidx',[],'mffkeysbackup',[]), event_new(rev_idx_opt(j):end)];
    end
end

openvar('event_new')
confirmation = input('Check empty lines are inserted at PR_. Continue? (y/n) ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('Continuing...');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%% check in variable event_new if empty lines are inserted in correct location (should
%%%%%%%%%%%%%% be instead of response (PR_) this time)

%introduce trigger names into struct
if length(event_new)==length(lst_xlsx_def)
    for i = 1:numel(event_new)
        event_new(i).type = lst_xlsx_def{i};
    end
else
    disp('Error: lengths of EGI triggers and excel triggers are unequal');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find out new position of empty rows
containsNumbers2 = cellfun(@(x) length(x)>0, {event_new.latency});

rowsToDelete = find(~containsNumbers2);

keepRows = true(size(event_new));
keepRows(rowsToDelete) = false;

% Delete empty rows
event_new = event_new(keepRows);

EEG.event = event_new;

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_trggrRepair.set'], 'filepath', EEG.filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 5: Epoch the data
% Before the prior and the update
timeBefore = 0.5;
% after the prior and the update
timeAfter = 1;
eventLabels={'PR_correct','PR_incorrect'};
updateEEG = pop_epoch(EEG,eventLabels,[-timeBefore timeAfter]); %this function looks at the type, not lable or code
EEG = updateEEG;
% Add epoch column to combined_files table and input epochs
numEpochs = EEG.trials;
epoch_numbers = (1:numEpochs)';
EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_epoched.set'], 'filepath', EEG.filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 6: Independent component analysis
% save the output of the previous step as a different
% variable - to be used in the next step

% If you are about to run a fresh ICA, but there are manually
% selected ICs from a previous iteration in this directory,
% then delete that manual ICs file
icFilefndr = dir([PROC_GROUP 'ICARemoved_*']);
if ~isempty(icFilefndr)
    delete([PROC_RS icFilefndr(1).name]);
end
% Load data

% Assuming EEG.data is a 3D array
[channels, timepoints, epochs] = size(EEG.data);
% Reshape the data into a 2D matrix, where each column represents time points across all epochs
EEGdata2D = reshape(EEG.data, channels, timepoints * epochs);
% Calculate the rank of the data
% should be roughly 128 - 1 minus # of channels removed
threshold = .04;
ChannelValues = svd(EEGdata2D); %svd(EEG.data);
inputRank = sum (ChannelValues(:) > threshold);
% run ica
EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',inputRank)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ICA.set'], 'filepath', EEG.filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%   Automatic ICA rejection

% assumes EEGLAB is preloaded on this computer
% [ALLEEG,~,~] = eeglab('rebuild');

EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, [0 0; ...     % Brain
    0.5 1; ...   % Muscle
    0.5 1; ...   % Eye
    0.5 1; ...   % Heart
    0.8 1; ...     % Line Noise
    0.5 1; ...     % Channel Noise
    0 0]);     % Other

% Identify rejected ICs
rejected_ICs = find(EEG.reject.gcompreject > 0);

% Remove the ICs flagged by ICLabel
EEG = pop_subcomp(EEG, [], 0);

%%%%%%%%%%%%%%  ONLY if you want to do a manual ICA rejection too

doManualICARejection = input('Do you want to perform manual ICA rejection? (yes/no): ', 's');

if strcmpi(doManualICARejection, 'yes')

    pop_viewprops(EEG,0);

    %write down bad components
    bad_components_add = [8,10,11,13,18,21,61,64];

    EEG.bad_components_add = sort(bad_components_add);
    EEG = epi_log(@pop_subcomp, EEG, EEG.bad_components_add, 1); % remove the components
    EEG = eeg_checkset(EEG, 'ica'); % update EEG.icaact

else

    disp('Skipping manual ICA rejection.');

end

%%%%%%%%%%%%%%  OTHERWISE: continue below

% Save removed ICAs
rejected_ICs_name = [EEG.filepath '\' 'ICS_rejected.mat'];
save(rejected_ICs_name,"rejected_ICs");
EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ICArej.set'], 'filepath', EEG.filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 9: Baseline correction
baselinePeriod = [-500 -300];

% Identify baseline time points
baselineTimePoints = EEG.times >= baselinePeriod(1) & EEG.times <= baselinePeriod(2);

% Calculate mean baseline for each channel
baselineValues = mean(EEG.data(:, baselineTimePoints, :), 2);

% Subtract baseline from each data point
EEG.data = EEG.data - repmat(baselineValues, 1, size(EEG.data, 2), 1);

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_baselineCorr.set'], 'filepath', EEG.filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 10: reject epoched data

pop_eegplot(EEG, 1, 1, 1); %have a look and check the triggers

confirmation = input('Do you want to close the trigger window? (y/n): ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('Window was closed closed.');
else
    error('Processing discontinued.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert struct array to table for easier manipulation
dataTable = struct2table(EEG.event);

% Group indices based on 'Age' column
[~, ~, groupIndices] = unique(dataTable.epoch);

% Create a cell array to store indices of each group
indicesByAge = accumarray(groupIndices, (1:numel(dataTable.epoch))', [], @(x) {x});

lst_stim = {EEG.event.type};
epochsToReject = [];

for i = 1:length(indicesByAge)
    lst_sm=lst_stim(indicesByAge{i});
    pr = cellfun(@(x) contains(x,'PR_'), lst_sm);
    idx_pr = find(pr);
    nm_pr = sum(pr);
    max_idx = max(idx_pr);
    if nm_pr>1
        % Define the indices of epochs to reject
        epochsToReject = [epochsToReject,i]; % Example indices of epochs to reject
    end
end

%%%%%%%%%%%%%% ONLY run below code if variable is NOT empty

if ~isempty(epochsToReject)
    % Reject epochs
    EEG = pop_rejepoch(EEG, epochsToReject, 0);

    % Save the dataset
    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_epochRej.set'], 'filepath', EEG.filepath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 11: define trials, correct vs incorrect
pop_eegplot(EEG, 1, 1, 1);
%load trial correctness information

answ = {EEG.event.type};
idx_ep_cor = [];
idx_ep_icor = [];
count=1; %this is the count for the number of epochs
for i = 1:numel(answ)
    if all(contains(answ{i},'PR_correct'))
        i;
        count;
        idx_ep_cor = [idx_ep_cor,count];
        count=count+1;
        answ{i};
    elseif all(contains(answ{i},'PR_incorrect'))
        i;
        count;
        idx_ep_icor = [idx_ep_icor,count];
        count=count+1;
        answ{i};
    end
end

EEG_corr = EEG.data([129],:,idx_ep_cor); %6:FCz, 11:Fz, 129:Cz
EEG_icorr = EEG.data([129],:,idx_ep_icor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
fs = 200;  % Sampling frequency in Hz

% Plot parameters
timeAxis = ((0:size(EEG_icorr,2)-1) / fs)-0.5;  % Time axis in seconds
% Plot butterfly plot for all electrodes
figure;
plot(timeAxis, squeeze(mean(EEG_icorr, 3)),'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Event-Related Potentials (ERPs) incorrect trials');
%legend(cellstr(num2str((1:size(EEG_icorr,1))')));  % Add electrode labels to legend
grid on;
filename_incorr = [EEG.filepath '\' EEG.filename(1:end-4) '_incorrect' , '.png'];
saveas(gcf, filename_incorr, 'png');

% Plot parameters
timeAxis = ((0:size(EEG_corr,2)-1) / fs)-0.5;  % Time axis in seconds
% Plot butterfly plot for all electrodes
figure;
plot(timeAxis, squeeze(mean(EEG_corr, 3)),'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Event-Related Potentials (ERPs) correct trials');
%legend(cellstr(num2str((1:size(EEG_icorr,1))')));  % Add electrode labels to legend
grid on;
filename_corr = [EEG.filepath '\' EEG.filename(1:end-4) '_correct' , '.png'];
saveas(gcf, filename_corr, 'png');

%plot difference error-correct

figure;
plot(timeAxis, squeeze(mean(EEG_icorr, 3))-squeeze(mean(EEG_corr, 3)),'r', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Event-Related Potentials (ERPs) correct, incorrect and difference');
%legend(cellstr(num2str((1:size(EEG_icorr,1))')));  % Add electrode labels to legend
grid on;
hold on;
plot(timeAxis, squeeze(mean(EEG_corr, 3)),'b', 'LineWidth', 0.5);
plot(timeAxis, squeeze(mean(EEG_icorr, 3)),'k', 'LineWidth', 0.5);
% Create text-only legend with associated colors
legendEntries = {'Difference', 'Correct', 'Incorrect'};
legendColors = {'r', 'b', 'k'};
for i = 1:numel(legendEntries)
    text(1, i, legendEntries{i}, 'Color', legendColors{i});
end
hold off;
set(gca, 'YDir', 'reverse');
filename_ERN = [EEG.filepath '\' EEG.filename(1:end-4) '_ERN' , '.png'];
saveas(gcf, filename_ERN, 'png');

confirmation = input('Do you want to close all figures? (y/n): ', 's');

% Check the user's response
if strcmpi(confirmation, 'y')
    close all;
    disp('All figures are closed.');
else
    error('Processing discontinued.');
end

EEG.corrIDX = EEG_corr
EEG.icorrIDX = EEG_icorr

%ERN calculated as mean value over 3 electrodes (Cz,Fz,FCz) and
%between 0 and 100ms (results in datapoints 100-120)
ERN_all = mean(mean(EEG_icorr, 3),1)
ERN = mean(ERN_all(100:120))
EEG.ERN = ERN

CRN_all = mean(mean(EEG_corr, 3),1)
CRN =  mean(CRN_all(100:120))
EEG.CRN = CRN

dERN = ERN-CRN
EEG.dERN = dERN

EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ERN_vals.set'], 'filepath', EEG.filepath);
