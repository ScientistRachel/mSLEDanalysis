% This script find dF/F for calcium fluorescence in an ROI around the SLED
% region.  Peak, neighbor, and persistence values are plotted & saved.
% dF/F is calculated with respect to the first frame of the time lapse.
%
% This script assumes the user has previously run find_tip_all.m
%
% Outputs: 
% (1) Figure showing wound ROI on top of original data
% (2) Plot showing peak, neighbor, and persistence dF/F
% (3) *_dF_F.mat file containing dF/F data for later plotting

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

% ROI region
dist_wound = 50; % Ignore non-cell regions that were not caused by the tip by ignoring blank space any further than this (in pixels) from the tip.
% Persistance definition
time_meas = 6; %in minutes, where to measure persistence

% Options
over_write = 1; % Overwrite files if a .mat file for this image already exists? 0 = no, 1 = yes
% Plotting Options
peak_lim = []; % Optional parameter to set the y-axis on the peak dF/F graph
pers_lim = []; % Optional parameter to set the y-axis on the persistence graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN ANALYSIS

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir filesep 'TimeAnalysisFigures' filesep 'IndidivualExperiments' filesep];
if ~exist([savedir filesep],'file')
    mkdir(savedir)
end

for mm = 1:length(directories)
    directory = directories{mm};
    
    %%%%% 2018/11/19 Add save tag based on directory so duplicate control names don't get confused
    savetag = strsplit(directory,'\');
    savetag = savetag{end-1};
 
    list = dir([directory '*.tif']);

    for jj = 1:length(list)
        imname = list(jj).name(1:end-4);
        
        if ~exist([savedir filesep savetag '_' imname '_dF_F.mat'],'file') || over_write

            load([mainDir imname '_tipFind.mat'])

            image_info = imfinfo([directory imname '.tif'],'tif');
            N_im = length(image_info)/2;

            % Boundaries of the ROI
            y_max = round(median(y_tip(~isnan(y_tip))))+dist_wound;
            y_min = round(median(y_tip(~isnan(y_tip))))-dist_wound;

            % Preallocate storage
            wound_int = NaN*ones(N_im,1);
            neigh_int = wound_int;
            % Loop through each frame of the time lapse
            for curr_im = 1:N_im

                im_Ca = imread([directory imname '.tif'],'tif',curr_im);

                % Only allow the wound to be near the tip
                im_w = im_Ca;
                im_w(1:y_min,:) = NaN;
                im_w(y_max:end,:) = NaN;

                % Neighbors = everything not near the wound
                im_n = im_Ca;
                im_n(y_min+1:y_max-1,:) = NaN;

                neigh_int(curr_im) = sum(sum(double(im_n)));
                wound_int(curr_im) = sum(sum(double(im_w)));

            end

            % Load in the original brightfield image for comparison to ROI
            im_bright = imread([directory imname '.tif'],'tif',curr_im+N_im);
            
            figure(1)
            imshow([imadjust(im_Ca) imadjust(im_bright)])
            hold on
            plot([1 size(im_Ca,2)*2 size(im_Ca,2)*2 1],[y_min y_min y_max y_max],'LineWidth',2)
            hold off
            saveas(gcf,[savedir filesep savetag '_' imname '_WoundROI.png'],'png')

            figure(2)
            set(gcf,'Position',[10 300 1120 420]) % This makes the figure larger, but might not be appropriate for all screens.
            time_vec = (1:N_im)*t_scale; %Time in seconds

            subplot(2,3,[1 2 4 5])
            plot(time_vec,(wound_int-wound_int(1))/wound_int(1),'Color','k','LineWidth',2)
            hold on
            plot(time_vec,(neigh_int-neigh_int(1))/neigh_int(1),'Color',0.6*[1 1 1],'LineWidth',2)
            plot(time_vec,0*time_vec,'--','Color',0.2*[1 1 1])
            hold off
            xlabel('Time (s)','FontSize',16)
            ylabel('\DeltaF/F','FontSize',16)
            box off
            legend('Edge','Neighbors')
            legend boxoff
            set(gca,'FontSize',16)

            peak_wound = max((wound_int-wound_int(1))/wound_int(1));
            peak_neigh = max((neigh_int-neigh_int(1))/neigh_int(1));

            subplot(2,3,3)
            b = bar(1,peak_wound);
            hold on
            b2 = bar(2,peak_neigh);
            hold off
            set(b,'FaceColor',[0 0 0],'EdgeColor','none')
            set(b2,'FaceColor',0.6*[1 1 1],'EdgeColor','none')
            set(gca,'XTick',1:2,'XTickLabel',{'Edge','Neighbors'},'FontSize',16,'XTickLabelRotation',90)
            ylabel('Peak \DeltaF/F','FontSize',16)
            title('Response')
            xlim([0.5 2.5])        
            if ~isempty(peak_lim)
                ylim(peak_lim)
            end


            % Find the persistence
            t400 = find(time_vec/60 == time_meas);
            if ~isempty(t400) % If there aren't enough frames, skip this
                wound400 = (wound_int(t400)-wound_int(1))/wound_int(1);
                neigh400 = (neigh_int(t400)-neigh_int(1))/neigh_int(1);

                subplot(2,3,6)
                b = bar(1,wound400);
                hold on
                b2 = bar(2,neigh400);
                hold off
                set(b,'FaceColor',[0 0 0],'EdgeColor','none')
                set(b2,'FaceColor',0.6*[1 1 1],'EdgeColor','none')
                set(gca,'XTick',1:2,'XTickLabel',{'Edge','Neighbors'},'FontSize',16,'XTickLabelRotation',90)
                ylabel(['DeltaF/F at ' num2str(time_meas) ' min'],'FontSize',16)
                title('Persistence')
                xlim([0.5 2.5])
                if ~isempty(pers_lim)
                    ylim(pers_lim)
                end
            end

            % Save the summary figure
            saveas(gcf,[savedir filesep savetag '_' imname '_dF_F.png'],'png')

            % Save dF/F data for later use
            curr_label = data_names{mm,2}{jj};
            save([savedir filesep savetag '_' imname '_dF_F.mat'],'wound_int','neigh_int',...
                'r_scale','t_scale','imname','N_im','time_vec','dist_wound',...
                'y_max','y_min','peak_wound','peak_neigh','wound400','neigh400',...
                'curr_label')

            clear wound_int neigh_int
        
        end

    end

end

close all
disp('Batch Complete')

