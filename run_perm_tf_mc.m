addpath(genpath( '/work/users/m/a/magcam/GAD/scripts/cluster_perm/plot_perm/AH_toolbox/'));
addpath(genpath( '/work/users/m/a/magcam/GAD/scripts/cluster_perm/plot_perm/Permutation_scripts/'));

dat_cor = load('/work/users/m/a/magcam/GAD/cluster_input/TF_correct_30ptps.mat');
dat_cor_3d = dat_cor.mat_end;

dat_icor = load('/work/users/m/a/magcam/GAD/cluster_input/TF_incorrect_30ptps.mat');
dat_icor_3d = dat_icor.mat_end;


thisCondDataT(:,1,:,:) = dat_cor_3d(:,:,:);
thisCondDataT(:,2,:,:) = dat_icor_3d(:,:,:);

ROI = 'rPC'

fdsRatio = 1;
tdsRatio = 1;

numIterations = 1000; % normally 1000, at sig level 0.05, 500 is the fewest number that would be useful
sigOptions = struct('onlyPos',0,'thresholdType','size');
minClusterSize = 20; 

permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); % Alpha threshol variable is the first significant threshold, for defining significance of individual values into the cluster

% 5. Call on AH_plotPerm to start permutation analysis
tic
[analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondDataT,{1,2},permutationOptions);
toc % readout how long perm test takes

perm.ROI = analysisStruct;
permMask.ROI = sigMaskInterp;
savePermName = append([ ROI '_perm_StimOnset.mat']);
savePermMaskName = append([ROI '_permMask_StimOnset.mat']);

% -- save variables 
save(append(['/work/users/m/a/magcam/GAD/cluster_output/tf2' savePermName]), 'perm');
save(append(['/work/users/m/a/magcam/GAD/cluster_output/tf2' savePermMaskName]), 'permMask');


lowerBoundFreq = 2;
upperBoundFreq = 58;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
FREQ = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);


% Assuming FREQ is your non-linearly spaced frequency vector
figure;
imagesc(-100/200:1/200:200/200, 1:length(FREQ), perm.ROI.real.p');
set(gca, 'YDir', 'normal');

% Set y-ticks to match frequency values
yticks = round(linspace(1, length(FREQ), 10)); % Pick 10 evenly spaced y-tick positions
yticklabels = arrayfun(@(x) sprintf('%.1f', FREQ(x)), yticks, 'UniformOutput', false); % Labels from FREQ
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

title('Time-Frequency Spectrogram p-values cluster permutation');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colormap(flipud(brewermap([], 'RdBu')));
colorbar;
caxis([0 0.05])



