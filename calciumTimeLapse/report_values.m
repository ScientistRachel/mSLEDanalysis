% This script combines dF/F information vs time and distance and displays
% summary statistics in the command windows
%
% This script assumes the user has previously find_dF_all_v2.m and
% find_dF_dist_all_v2.m
%
% Output: Summary statistics on the command line

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\';

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

%% Display Peak and Persistence Information

WoundPeak = NaN*ones(length(labels),N_total);
NeighPeak = WoundPeak;
Wound400 = WoundPeak; Neigh400 = WoundPeak;

list = dir([mainDir 'TimeAnalysisFigures' filesep 'IndidivualExperiments' filesep '*dF_F.mat']);

for kk = 1:length(list)
   
    load([mainDir 'TimeAnalysisFigures' filesep 'IndidivualExperiments' filesep list(kk).name])

    curr_label = find(strcmp(labels,curr_label));
    
    WoundPeak(curr_label,kk) = peak_wound;
    NeighPeak(curr_label,kk) = peak_neigh;
    Wound400(curr_label,kk) = wound400;
    Neigh400(curr_label,kk) = neigh400;
    
end

WoundMean = nanmean(WoundPeak,2);
WoundStd = nanstd(WoundPeak,[],2);

WoundMean400 = nanmean(Wound400,2);
WoundStd400 = nanstd(Wound400,[],2);

NeighMean = nanmean(NeighPeak,2);
NeighStd = nanstd(NeighPeak,[],2);

NeighMean400 = nanmean(Neigh400,2);
NeighStd400 = nanstd(Neigh400,[],2);

disp(' ')
disp('Summary of Peak Data:')
disp(' ')

T_peak = table(WoundMean,WoundStd,NeighMean,NeighStd,...
    'RowNames',labels,'VariableNames',{'MeanWound' 'SDWound' 'MeanNeighbors' 'SDNeighbors'});

disp(T_peak)
    
disp(' ')
disp('Summary of Persistence Data (dF/F at 6 min):')
disp(' ')

T_peak = table(WoundMean400,WoundStd400,NeighMean400,NeighStd400,...
    'RowNames',labels,'VariableNames',{'MeanWound' 'SDWound' 'MeanNeighbors' 'SDNeighbors'});

disp(T_peak)


disp(' ')
disp('N values')
disp(' ')

N = sum(~isnan(WoundPeak'));

T_N = table(N','RowNames',labels,'VariableNames',{'N'});
disp(T_N)


%% Display Distance Information

peakDist = NaN*ones(length(labels),N_total,3);

list = dir([mainDir 'DistAnalysisFigures' filesep 'IndidivualExperiments' filesep '*dF_F_dist.mat']);

% Gather the data together
for kk = 1:length(list)
   
    load([mainDir 'DistAnalysisFigures' filesep 'IndidivualExperiments' filesep list(kk).name])

    curr_label = find(strcmp(labels,curr_label));
    
    peakDist(curr_label,kk,:) = peak_dist;
    
end

Dist1 = peakDist(:,:,1)';
Dist2 = peakDist(:,:,2)';
Dist3 = peakDist(:,:,3)';

Dist150_Mean = nanmean(Dist1);
Dist150_Std = nanstd(Dist1);

Dist250_Mean = nanmean(Dist2);
Dist250_Std = nanstd(Dist2);

Dist500_Mean = nanmean(Dist3);
Dist500_Std = nanstd(Dist3);

disp('Summary of Distance Data:')
disp(' ')

T_dist = table(Dist150_Mean(:),Dist150_Std(:),Dist250_Mean(:),Dist250_Std(:),Dist500_Mean(:),Dist500_Std(:),...
    'RowNames',labels,'VariableNames',{'Mean150' 'SD150' 'Mean250' 'SD250' 'Mean500' 'SD500'});

disp(T_dist)

