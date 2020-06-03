% This script combines dF/F information from multiple replicates into
% figures with standard deviation error bars
%
% This script assumes the user has previously find_dF_dist_all_v2.m
%
% Outputs: Figures showing dF/F as a function of distance for each 
% experiment type
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
ylim_user = [-1 8]; %used on combined dF/F plots, set to [] to ignore
labels_use = {'MCF10A','PTEN^{-/-}'}; % If left empty, defaults to the labels used in do_directory_setup, but can be used to make 'prettier' labels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir filesep 'DistAnalysisFigures' filesep];
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
xAxis = cell(length(labels),length(directories));
maxDist = xAxis;

peakDist = NaN*ones(length(labels),length(directories),3);
threshDist = NaN*ones(length(labels),length(directories));

list = dir([mainDir 'DistAnalysisFigures\IndidivualExperiments' filesep '*dF_F_dist.mat']);

% Gather the data together
Lin = zeros(length(labels),1);
for kk = 1:length(list)
   
    load([mainDir 'DistAnalysisFigures\IndidivualExperiments' filesep list(kk).name])

    curr_label = find(strcmp(labels,curr_label));
    Lin(curr_label) = Lin(curr_label)+1;
    Lkk = Lin(curr_label);
    
    xAxis{curr_label,Lkk} = max_dist_vec;
    maxDist{curr_label,Lkk} = max_dist;
        
end

% Preallocate storage
allMeans = cell(length(labels),1);
allStd = allMeans;
allX = allMeans;

% Take averages of distance curves
for kk = 1:length(labels)
    
    slice = xAxis(kk,:);
    slice2 = maxDist(kk,:);
    
    % Find the shortest experiment to use for taking averages
    Lmin = cellfun('length',slice);
    Lmin = min(Lmin(:));
    
    count = 1;
    sliceX = NaN*ones(1,Lmin);
    sliceM = sliceX;
    for jj = 1:length(slice)
        if ~isempty(slice{jj})
            sliceX(count,:) = slice{jj}(1:Lmin);
            sliceM(count,:) = slice2{jj}(1:Lmin);
            count = count+1;
        end
    end
    
    if sum(sum(diff(sliceX))) > 0 && size(sliceX,1) > 1 % Sanity check
        error('Distances don''t add up')
    end
    
    sliceMean = mean(sliceM);
    sliceStd = std(sliceM);
    sliceN = sum(~isnan(sliceM));
    xMean = mean(sliceX);
    if size(sliceX,1) == 1 % If there wasn't anything to average over, keep the full vector
        xMean = sliceX;
        sliceMean = sliceM;
        sliceStd = NaN*sliceM;
        sliceN = sliceM./sliceM;
    end
    
    allMeans{kk} = sliceMean;
    allStd{kk} = sliceStd;
    allX{kk} = xMean;
    
    t_star = tinv(1-0.05/2,sliceN-1); % Not currently used -- could be used to turn the standard deviation error bars into 95% CI
    
    figure(1)
    h = boundedline(xMean,sliceMean,sliceStd,.../sqrt(sliceN).*t_star,...
        'cmap',[0 0 0 ; 0.6 0.6 0.6]);
    hold on
    plot(xMean,0*xMean,'--','Color',0.2*[1 1 1])
    hold off
    set(h,'LineWidth',2)
    xlabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
    ylabel('Peak \DeltaF/F','FontSize',16)
    box off
    set(gca,'FontSize',16)
    xlim([min(xMean) max(xMean)])
    title(labels_use{kk})
    if ~isempty(ylim_user)
        ylim(ylim_user)
    end
    
    saveas(gcf,[savedir labels{kk} '_dF_F_dist.tif'],'tif')
    
    close all
     
end  
    
