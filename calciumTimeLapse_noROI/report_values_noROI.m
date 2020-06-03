% This script combines dF/F information vs time and displays
% summary statistics in the command windows
%
% This script assumes the user has previously find_dF_all_noROI.m
%
% Output: Summary statistics on the command line

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

N_total = 0;
for kk = 1:size(data_names,1)
    N_total = N_total + length(data_names{kk,2});
end

%% Peak and Persistence

WoundPeak = NaN*ones(length(labels),N_total);
Wound400 = WoundPeak;

list = dir([mainDir 'TimeAnalysisFigures\IndidivualExperiments' filesep '*dF_F.mat']);

for kk = 1:length(list)
   
    load([mainDir 'TimeAnalysisFigures\IndidivualExperiments' filesep list(kk).name])

    curr_label = find(strcmp(labels,curr_label));
    
    WoundPeak(curr_label,kk) = peak_wound;
    Wound400(curr_label,kk) = wound400;
    
end

WoundMean = nanmean(WoundPeak,2);
WoundStd = nanstd(WoundPeak,[],2);

WoundMean400 = nanmean(Wound400,2);
WoundStd400 = nanstd(Wound400,[],2);

disp(' ')
disp('Summary of Peak Data:')
disp(' ')

T_peak = table(WoundMean,WoundStd,...
    'RowNames',labels,'VariableNames',{'MeanWound' 'SDWound'});

disp(T_peak)
    
disp(' ')
disp('Summary of Persistence Data (dF/F at 6 min):')
disp(' ')


T_peak = table(WoundMean400,WoundStd400,...
    'RowNames',labels,'VariableNames',{'MeanWound' 'SDWound'});

disp(T_peak)

disp(' ')
disp('N values')
disp(' ')

N = sum(~isnan(WoundPeak'));

T_N = table(N','RowNames',labels,'VariableNames',{'N'});
disp(T_N)