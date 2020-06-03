% This script find dF/F for calcium fluorescence in an image.
% Peak and persistence values are plotted & saved.
% dF/F is calculated with respect to the first frame of the time lapse.
%
% Outputs: 
% (1) Plot showing peak, neighbor, and persistence dF/F
% (2) *_dF_F.mat file containing dF/F data for later plotting

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

% Persistance definition
time_meas = 6; %in minutes, where to measure persistence

% Options
over_write = 1; % Overwrite files if a .mat file for this image already exists? 0 = no, 1 = yes
% Plotting Options
peak_lim = []; % Optional parameter to set the y-axis on the peak dF/F graph
pers_lim = []; % Optional parameter to set the y-axis on the persistence graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Code

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir filesep 'TimeAnalysisFigures' filesep 'IndidivualExperiments' filesep];
if ~exist([savedir filesep],'file')
    mkdir(savedir)
end

savetag = [];
for mm = 1:length(directories)
    directory = directories{mm};
    
    %%%%% 2018/11/19 Add saving so duplicate control names don't get confused
    savetag = strsplit(directory,'\');
    savetag = savetag{end-1};
 
    list = dir([directory '*.tif']);

    for jj = 1:length(list)
        imname = list(jj).name(1:end-4);
        
        if ~exist([savedir filesep savetag '_' imname '_dF_F.mat'],'file') || over_write

        image_info = imfinfo([directory imname '.tif'],'tif');
        N_im = length(image_info)/2;
        

        wound_int = NaN*ones(N_im,1);
        for curr_im = 1:N_im

            im_Ca = imread([directory imname '.tif'],'tif',curr_im);
            wound_int(curr_im) = sum(sum(double(im_Ca)));

        end
        
        time_vec = (1:N_im)*t_scale; %Time in seconds
                
        figure(2)
        set(gcf,'Position',[10 300 1120 420])

        subplot(2,3,[1 2 4 5])
        plot(time_vec,(wound_int-wound_int(1))/wound_int(1),'Color','k','LineWidth',2)
        hold on
        plot(time_vec,0*time_vec,'--','Color',0.2*[1 1 1])
        hold off
        xlabel('Time (s)','FontSize',16)
        ylabel('\DeltaF/F','FontSize',16)
        box off
        set(gca,'FontSize',16)
        
        peak_wound = max((wound_int-wound_int(1))/wound_int(1));
        
        subplot(2,3,3)
        b = bar(1,peak_wound);

        set(b,'FaceColor',[0 0 0],'EdgeColor','none')
        set(gca,'XTick',[],'FontSize',16,'XTickLabelRotation',90)
        ylabel('Peak \DeltaF/F','FontSize',16)
        title('Response')
        xlim([0 2])        
        if ~isempty(peak_lim)
            ylim(peak_lim)
        end
        
        
        t400 = find(time_vec/60 == time_meas);
        if ~isempty(t400)
            wound400 = (wound_int(t400)-wound_int(1))/wound_int(1);
            
            subplot(2,3,6)
            b = bar(1,wound400);

            set(b,'FaceColor',[0 0 0],'EdgeColor','none')
            set(gca,'XTick',[],'FontSize',16,'XTickLabelRotation',90)
            ylabel(['DeltaF/F at ' num2str(time_meas) ' min'],'FontSize',16)
            title('Persistence')
            xlim([0 2])
            if ~isempty(pers_lim)
                ylim(pers_lim)
            end
        end

        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F.png'],'png')
        
        curr_label = data_names{mm,2}{jj};

        save([savedir filesep savetag '_' imname '_dF_F.mat'],'wound_int',...
            'r_scale','t_scale','imname','N_im','time_vec',...
            'peak_wound','wound400',...
            'curr_label')
        
        clear wound_int
        
        end

    end

end

close all
disp('Batch Complete')