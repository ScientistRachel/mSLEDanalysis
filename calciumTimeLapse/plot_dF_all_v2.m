% This script combines dF/F information from multiple replicates into
% figures with standard deviation error bars
%
% This script assumes the user has previously find_dF_all_v2.m
%
% Outputs: Figure showing dF/F as a function of time for each experiment type
%
% Dependencies:
% Kelly Kearney (2020). boundedline.m (https://www.github.com/kakearney/boundedline-pkg), GitHub. Retrieved June 2, 2020.

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

% Plotting Options
ylim_user = [-.5 4]; %used on combined dF/F plots, set to [] to ignore
labels_use = {'MCF10A','PTEN^{-/-}'}; % If left empty, defaults to the labels used in do_directory_setup, but can be used to make 'prettier' labels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run analysis

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir filesep 'TimeAnalysisFigures' filesep];
if ~exist([savedir filesep],'file')
    mkdir(savedir)
end

if isempty(labels_use)
    labels_use = labels;
end
if length(labels_use) ~= length(labels)
    error('Bad Labeling')
end

% Preallocate storage for total data
WoundInt = cell(length(labels),1);
NeighInt = WoundInt;

list = dir([mainDir 'TimeAnalysisFigures' filesep 'IndidivualExperiments' filesep '*dF_F.mat']);

% Gather the data together
Lin = zeros(length(labels),1);
for kk = 1:length(list)
   
    load([mainDir 'TimeAnalysisFigures\IndidivualExperiments' filesep list(kk).name])

    curr_label = find(strcmp(labels,curr_label));
    Lin(curr_label) = Lin(curr_label)+1;
    Lkk = Lin(curr_label);
    
    WoundInt{curr_label}(:,Lkk) = (wound_int-wound_int(1))/wound_int(1);
    NeighInt{curr_label}(:,Lkk) = (neigh_int-neigh_int(1))/neigh_int(1);
    
end

if exist('time_vec','var')
    disp('time_vec exists')
else
    error('find the time_vec')
end

clear wound_int neigh_int peak_wound peak_neigh wound400 neigh400 dist_wound
clear N_im curr_label dF_file directory imname y_max y_min

% Take averages of time curves
WoundMean = NaN*ones(200,length(labels));
WoundStd = WoundMean; WoundN = WoundMean;

NeighMean = WoundMean;
NeighStd = WoundMean; NeighN = WoundMean;

for kk = 1:length(labels)
    
    IntNow = WoundInt{kk};
    IntNow = cutNaNdir(IntNow')';
    
    WoundMean(:,kk) = mean(IntNow,2);
    WoundStd(:,kk) = std(IntNow,[],2);
    WoundN(:,kk) = sum(~isnan(IntNow),2);
    
    IntNow = NeighInt{kk};
    IntNow = cutNaNdir(IntNow')';
    
    NeighMean(:,kk) = mean(IntNow,2);
    NeighStd(:,kk) = std(IntNow,[],2);
    
    t_star = tinv(1-0.05/2,WoundN(:,kk)-1); % Not currently used -- could be used to turn the standard deviation error bars into 95% CI
    
    figure(1)
    plot(time_vec,WoundMean(:,kk),'Color','k','LineWidth',2)
    hold on
    plot(time_vec,NeighMean(:,kk),'Color',0.6*[1 1 1],'LineWidth',2)
    plot(time_vec,0*time_vec,'--','Color',0.2*[1 1 1])
    h = boundedline(time_vec,WoundMean(:,kk),WoundStd(:,kk),..../sqrt(WoundN(:,kk)).*t_star,...
        time_vec,NeighMean(:,kk),NeighStd(:,kk),..../sqrt(WoundN(:,kk)).*t_star,...
        'cmap',[0 0 0 ; 0.6 0.6 0.6]);
    set(h,'LineWidth',2)
    hold off
    xlabel('Time (s)','FontSize',18)
    ylabel('\DeltaF/F','FontSize',18)
    box off
    legend('Edge','Neighbors')
    legend boxoff
    set(gca,'FontSize',18)
    title(labels_use{kk})
    xlim([min(time_vec(:)) max(time_vec(:))])
    if ~isempty(ylim_user)
        ylim(ylim_user)
    end
    
    saveas(gcf,[savedir labels{kk} '_dF_F_ErrorBarStd.tif'],'tif')

end

close all
