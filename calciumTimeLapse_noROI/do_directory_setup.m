% This file is the master list of what files contain what data for later
% figure plotting.  Saved as a .mat file to be called by other scripts.

clc, clear, close all

%%%%%% User Inputs: Location, Type, and Replicate Information for Each Experiment
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\ExampleOutput\'; % Location for output from this and other scripts

% Directory 1
data_names{1,1} = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\imagesTIF\'; % Directory containing images
data_names{1,2} = {'Add Control','Add H_2O_2','No Add Control'}; % Type of data (e.g. which cell line or which drug)
data_names{1,3} = repmat({'111418'},[1 3]); % Replicate information (which experiment, N = 1,2, etc. or date of experiment)

% Can add additional directories (2, 3, ...N) here:
% data_names{end+1,1} = 'add folder here';
% data_names{end,2} = 'add types here';
% data_names{end,3} = 'add replicates here';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Process inputs

% Automatically detect the names of different data types for later use
labels = data_names(:,2);
labels = [labels{:}];
labels = unique(labels);
%Check for bad data
bad_check = strcmp('BAD',labels);
labels(bad_check) = [];
clear bad_check

%Show the user a sanity check
disp('These are the types you have entered:')
for kk = 1:numel(labels), disp(labels{kk}), end
clear kk

if ~exist(mainDir,'file')
    mkdir(mainDir)
end

save([mainDir 'SLED_data_names.mat'])

disp(' ')