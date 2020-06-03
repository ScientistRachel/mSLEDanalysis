% Script analyzes the difference in PI staining between two images.
% Assumes that images are .tif files with paired file names such that there
% are files with the same name except for a designation of before and after
% (e.g. "Prescratch" and "Postscratch")

clc, clear, close all

%%%%% User Inputs

% Image directory and save directory
imDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\PI_staining\imagesTIF\PI\';
saveDir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\PI_staining\ExampleOutput\';

% Plot options
plot_on = 1; % When = 1, plots comparison between pre- and post-images.  When = 0, no plots are shown or saved
fontsize = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the analysis

if ~exist(saveDir,'file')
    mkdir(saveDir)
end

list = dir([imDir '*.tif']);
% Preallocate storage
N = NaN*ones(length(list)/2,1); % Variable storing summed pixel level change in intensities
for kk = 1:length(list)/2
    
    imPost = imread([imDir list(2*kk-1).name]);
    imPre = imread([imDir list(2*kk).name]);
    
    disp(['Comparing ' list(2*kk-1).name ' and ' list(2*kk).name])
    
    L = size(imPre,2);
    
    % Measure how many pixels are saturated (the lower, the better quality
    % this analysis will be)
    ReSca = [imPre imPost];
    disp(['     Percent Saturated Pixels: ' num2str(100*sum(sum(ReSca == 255))/L^2) '%'])
    
    imPre = ReSca(:,1:L);
    imPost = ReSca(:,L+1:end);

    diffIm = imPost-imPre;
    diffIm(diffIm < 0) = 0;
    n = sum(diffIm(:));
   
    if plot_on > 0 % Optionally, plot a comparison of the pre and post images
        figure(1)
        set(gcf,'Position',[10          50        1680         840]) % Can comment out or change this line changing figure size as approrpiate for your screen

        subplot(1,2,1)
        imshowpair(imPre,imPost)

        subplot(1,2,2)
        imshow(diffIm)

        title(['\DeltaI = ' num2str(n)])
        saveas(gcf,[saveDir 'ComparePrePost_NewPixels_' list(2*kk-1).name(1:end-19) '.png'],'png')
    end
    
    N(kk) = n;
    
end


%% Example Figures
% Labeling on these figures is based on the set of example images shared
% with this code and needs to be updated based on the actual data being
% analyzed.

% Figuring showing raw pixel value increase in PI staining
figure(2)
bar(N,'FaceColor',[0    0.6196    0.4510],'EdgeColor',0.5*[0    0.6196    0.4510],'LineWidth',2)
set(gca,'XTickLabel',{'0.2 kPa','8 kPa','64 kPa'},'FontSize',fontsize)
ylabel('\DeltaI (raw pixels)','FontSize',fontsize)
title('Increase in PI Staining')
box off
saveas(gcf,[saveDir 'RawPixelChangeComparison.png'],'png')

% Figure showing pixel values normalized to a control
figure(3)
bar(N/N(1),'FaceColor',[0    0.6196    0.4510],'EdgeColor',0.5*[0    0.6196    0.4510],'LineWidth',2)
hold on
plot([0 4],[1 1],'--','Color',0.6*[1 1 1])
hold off
set(gca,'XTickLabel',{'0.2 kPa','8 kPa','64 kPa'},'FontSize',fontsize)
ylabel('\DeltaI / \DeltaI(0.2 kPa)','FontSize',fontsize)
title('Increase in PI Staining')
box off
saveas(gcf,[saveDir 'NormalizedPixelChangeComparison.png'],'png')


%%
close all