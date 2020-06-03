% This script calculates at dF/F as a function of distance from the SLED
% pipette tip. dF/F is calculated with respect to the first frame of the
% time lapse.
%
% This script assumes the user has previously run find_tip_all.m
%
% Outputs: 
% (1) Kymographs of (a) raw intensity and (b) dF/F
% (2) Curves showing dF/F as a function of (a) time and (b) distance
% (3) *_dF_F_dist.mat file containing dF/F data for later plotting

clc, clear, close all

%%%%%% User Inputs
% Where should the data be saved?
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\';

% Image parameters
t_scale = 2; %seconds
r_scale = 2.485; %um/pixel

% Analysis parameter: which distance(s) should the max(dF/F) be calculated?
max_dist_meas = [ceil(150/r_scale) ceil(250/r_scale) ceil(500/r_scale)]; % 150, 250, and 500 um, converted to nearest pixel value

% Options
over_write = 1; % Overwrite files if a .mat file for this image already exists? 0 = no, 1 = yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN ANALYSIS

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

savedir = [mainDir 'DistAnalysisFigures' filesep 'IndidivualExperiments'];
if ~exist([savedir filesep],'file')
    mkdir(savedir)
end

for mm = 1:length(directories)
    directory = directories{mm};
    
    %%%%% 2018/11/19 Add saving so duplicate control names don't get confused
    savetag = strsplit(directory,'\');
    savetag = savetag{end-1}; 
    
    list = dir([directory '*.tif']);

    for jj = 1:length(list)
        imname = list(jj).name(1:end-4);
        
        if ~exist([savedir filesep savetag '_' imname '_dF_F_dist.mat'],'file') || over_write

        load([mainDir imname '_tipFind.mat'])

        image_info = imfinfo([directory imname '.tif'],'tif');
        N_im = length(image_info)/2;        

        % Average the information from each frame onto the y-axis
        yproj = NaN*ones(image_info(1).Height,N_im);
        for curr_im = 1:N_im
            im_Ca = imread([directory imname '.tif'],'tif',curr_im);            
            yproj(:,curr_im) = mean(double(im_Ca),2);
        end
        % dF/F calculated with respect to first frame
        % (normalization is dependent on y-axis information)
        yproj_dF = (yproj - repmat(yproj(:,1),[1 size(yproj,2)]))./repmat(yproj(:,1),[1 size(yproj,2)]);
        
        % Useful info for plotting
        time_vec = (1:N_im)*t_scale; %Time in seconds
        dist_vec = (1:size(im_Ca,2))*r_scale; % Distance in um
        centerline = round(median(y_tip(~isnan(y_tip)))); % Location of SLED pipette tip        
        dist_cent = dist_vec - (centerline*r_scale); % Distance from SLED pipette tip
        
        % Raw intensity kymograph
        figure(1)
        imagesc([min(time_vec) max(time_vec)],[min(dist_cent) max(dist_cent)],yproj)
        xlabel('Time (s)','FontSize',16)
        ylabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        set(gca,'FontSize',16)
        h = colorbar;
        set(get(h,'Label'),'String','Raw Fluorescence Intensity','FontSize',16)
        saveas(gcf,[savedir filesep savetag '_' imname '_IntensityKymograph.png'],'png')
        
        % dF/F kymograph
        figure(2)
        imagesc([min(time_vec) max(time_vec)],[min(dist_cent) max(dist_cent)],yproj_dF)
        xlabel('Time (s)','FontSize',16)
        ylabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        set(gca,'FontSize',16)
        h = colorbar;
        caxis manual
        caxis([0 4.5])
        set(get(h,'Label'),'String','\DeltaF/F','FontSize',16)
        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F_Kymograph.png'],'png')
        
        % dF/F curves colored by time
        figure(3)
        cmap_kymo = parula(size(yproj_dF,2));
        for kk = 1:size(yproj_dF,2)
            plot(dist_cent,yproj_dF(:,kk),'color',cmap_kymo(kk,:))
            if kk == 1
                hold on
            end
        end
        hold off
        ylabel('\DeltaF/F','FontSize',16)
        xlabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        set(gca,'FontSize',16)
        xlim([min(dist_cent) max(dist_cent)])
        h = colorbar;
        caxis manual
        caxis([min(time_vec) max(time_vec)])
        set(get(h,'Label'),'String','Time (s)','FontSize',16)
        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F_kymoCurves_colorTime.png'],'png')
        
        % dF/F curves colored by distance
        figure(4)
        n_dists = min(size(yproj_dF,1)-centerline,centerline);
        cmap_kymo = parula(n_dists);
        dist_vals = NaN*ones(n_dists,size(yproj_dF,2));
        for kk = 0:n_dists-1
            slice = yproj_dF([centerline-kk centerline+kk],:);
            plot(time_vec,mean(slice),'color',cmap_kymo(kk+1,:))
            if kk == 1
                hold on
            end
            dist_vals(kk+1,:) = mean(slice);
        end
        hold off
        ylabel('\DeltaF/F','FontSize',16)
        xlabel('Time (s)','FontSize',16)
        set(gca,'FontSize',16)
        xlim([min(time_vec) max(time_vec)])
        h = colorbar;
        caxis manual
        caxis([0 (n_dists-1)*r_scale])
        set(get(h,'Label'),'String',['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F_kymoCurves_colorDistance.png'],'png')

        % Calculate the maximum at each distance
        max_dist = max(dist_vals,[],2);
        max_dist_vec = (0:(n_dists-1))*r_scale;
        
        
        figure(5)
        set(gcf,'Position',[10 300 1000 500]) % Make the figure bigger; might not work for all monitors
        subplot(1,3,1:2) % Curve of max dF/F vs distance
        plot(max_dist_vec,max_dist,'Color','k','LineWidth',2)
        xlabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        ylabel('Peak \DeltaF/F','FontSize',16)
        box off
        set(gca,'FontSize',16)
        xlim([min(max_dist_vec) max(max_dist_vec)])
        
        subplot(1,3,3) % Max dF/F at user chosen distances
        try
            peak_dist = max_dist(max_dist_meas);
            x_tick_label = num2str(round(max_dist_vec(max_dist_meas)'));
        catch % If the above doesn't work, there isn't enough distance from the scratch for the last value
            error('Not enough distance from scratch')
        end
        b = bar(peak_dist);
        set(b,'FaceColor',[0 0 0],'EdgeColor','none')
        set(gca,'XTickLabel',x_tick_label,'FontSize',16)
        ylabel('Peak \DeltaF/F','FontSize',16)
        xlabel(['Distance from Scratch (' char(181) 'm)'],'FontSize',16)
        xlim([0 length(max_dist_meas)]+.5)

        saveas(gcf,[savedir filesep savetag '_' imname '_dF_F_dist.png'],'png')
        
        % Save data for later plotting
        curr_label = data_names{mm,2}{jj};
        save([savedir filesep savetag '_' imname '_dF_F_dist.mat'],...
            'r_scale','t_scale','imname','N_im','time_vec','dist_vec',...
            'centerline','dist_cent','peak_dist','yproj_dF',...
            'max_dist','max_dist_vec','max_dist_meas',...
            'curr_label')
            
        clear b cmap_kymo curr_im h im_Ca image_info slice test_thresh tip_mask x_tip y_tip
        
        end

    end

end

close all
disp('Batch Complete')