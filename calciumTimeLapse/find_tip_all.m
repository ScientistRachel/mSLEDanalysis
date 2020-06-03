% This script uses the brightfield information to find the pipette tip used
% in SLED stimulation.  The tip of the pipette will be used to define
% regions of interest (ROIs) in later analysis.
%
% This script assumes the user has previously run do_directory_setup.m
%
% Output: _tipFind.mat files with x & y locations of pipette over time

clc, clear, close all

%%%%%% User inputs
% Location for output from this and other scripts
mainDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse\ExampleOutput\'; 
% Options:
show_images = 1; % If set to 1, will plot the found tip location.  Recommend to set to 1 initially to confirm accuracy
overwrite = 0; % Overwrite files if a .mat file for this image already exists? 0 = no, 1 = yes
% Image analysis parameters:
% (Note these parameters are robust to multiple replicates and should be
% held constanst for compared experiments):
tip_thresh = 1000; % BW threshold for the tip darkness
strel_disk = 10; % Structuring element to clean up the image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN TIP FINDING

load([mainDir 'SLED_data_names.mat'])
directories = data_names(:,1);

for mm = 1:length(directories)
    directory = directories{mm};
    list = dir([directory '*.tif']);

    for jj = 1:length(list)
        imname = list(jj).name(1:end-4);
        disp(imname)
        
        if ~exist([mainDir imname '_tipFind.mat'],'file') || overwrite

            image_info = imfinfo([directory imname '.tif'],'tif');
            N_im = length(image_info)/2;

            x_tip = NaN*ones(N_im,1);
            y_tip = x_tip;
            tip_mask = NaN*ones(image_info(1).Height,image_info(1).Width,N_im);
            for kk = 1:N_im

                im_tip = imread([directory imname '.tif'],'tif',N_im+kk); % brightfield is second half of the channel
                bw_tip = im_tip < tip_thresh;
                bw_tip = imopen(bw_tip,strel('disk',strel_disk));
                [y,x] = find(bw_tip);
                tip = find(x == min(x(:)));

                if show_images > 0 && kk < 50 
                    imshow(imadjust(im_tip))
                    hold on
                    plot(mean(x(tip)),mean(y(tip)),'o','Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',10)
                    drawnow
                    hold off
                end

                x_tip(kk) = mean(x(tip));
                y_tip(kk) = mean(y(tip));
                
                tip_mask(:,:,kk) = bw_tip;

            end

            save([mainDir imname '_tipFind.mat'],'x_tip','y_tip','tip_mask')
        
        end

    end

end

disp(' ')
close all
disp('Batch Completed')