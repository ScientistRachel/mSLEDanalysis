% This script plots kymographs, similar to find_dF_dist_all_v2.m, but with
% larger fonts and prettier labels.
%
% This script assumes the user has previously run find_tip_all.m
%
% Outputs: Kymographs of dF/F

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

% Options
over_write = 1; % Overwrite files if a .png file for this experiment already exists? 0 = no, 1 = yes
% Plotting Options
fontsize = 20;
labels_use = {'MCF10A','PTEN^{-/-}'}; % If left empty, defaults to the labels used in do_directory_setup, but can be used to make 'prettier' labels.
caxisLim = [0 5]; % Consistant colorbar axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir 'Kymographs_LargeFont' filesep];
if ~exist([savedir filesep],'file')
    mkdir(savedir)
end

if isempty(labels_use)
    labels_use = labels;
end
if length(labels_use) ~= length(labels)
    error('Bad Labeling')
end


for mm = 1:length(directories)
    directory = directories{mm};
    
    %%%%% 2018/11/19 Add directory to saving so duplicate control names don't get confused
    savetag = strsplit(directory,'\');
    savetag = savetag{end-1};
    
    list = dir([directory '*.tif']);

    for jj = 1:length(list)
        imname = list(jj).name(1:end-4);
        label_now = labels_use(strcmp(labels,data_names{mm,2}{jj}));
        
        if ~exist([savedir filesep imname '_dF_F_Kymograph.png'],'file') || over_write

        load([mainDir imname '_tipFind.mat'])

        image_info = imfinfo([directory imname '.tif'],'tif');
        N_im = length(image_info)/2;        

        % Average over x-dimension to only leave y-dimension information
        yproj = NaN*ones(image_info(1).Height,N_im);
        for curr_im = 1:N_im
            im_Ca = imread([directory imname '.tif'],'tif',curr_im);            
            yproj(:,curr_im) = mean(double(im_Ca),2);
        end
        % dF/F normalized to first frame
        yproj_dF = (yproj - repmat(yproj(:,1),[1 size(yproj,2)]))./repmat(yproj(:,1),[1 size(yproj,2)]);
        
        % Useful for plots
        time_vec = (1:N_im)*t_scale; %Time in seconds
        dist_vec = (1:size(im_Ca,2))*r_scale; % Distance in um
        centerline = round(median(y_tip(~isnan(y_tip))));
        dist_cent = dist_vec - (centerline*r_scale);
        
        figure(1)
        imagesc([min(time_vec) max(time_vec)],[min(dist_cent) max(dist_cent)],yproj_dF)
        xlabel('Time (s)','FontSize',fontsize)
        ylabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',fontsize)
        
        %Clunky way to get reasonable y-tick
        max_y = max(abs(dist_cent));
        ticks = sort([0:200:max_y -200:-200:-max_y]);
        set(gca,'YTick',ticks)
        
        set(gca,'FontSize',fontsize)
        h = colorbar;
        caxis manual
        caxis(caxisLim)
        set(get(h,'Label'),'String','\DeltaF/F','FontSize',fontsize)
        
        title(label_now)
        
        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F_Kymograph.png'],'png')
        
        end

    end

end

close all
disp('Batch Complete')