% wipe out all previous data
clear all; clc; close all;

scripts = '/work/users/m/a/magcam/GAD/scripts/';
addpath(scripts)

root_raw = '/work/users/m/a/magcam/GAD/EEG_mats/';

% Define the top-level directory and the pattern
topLevelDir = root_raw; % Replace with your directory path
filePattern = '**/EEG_ptp*.mat'; % Replace with your file pattern (e.g., '*.txt' or '*.mat')

% Search for all matching files
files = dir(fullfile(topLevelDir, filePattern));
% Preallocate cell array to store full paths
filePaths = cell(length(files), 1);
% Store the full paths of the matching files in the cell array
for k = 1:length(files)
    filePaths{k} = fullfile(files(k).folder, files(k).name);
end


for i = 19:length(filePaths)
    PACz_all3_MC(filePaths{i})
end

