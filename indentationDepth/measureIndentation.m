% This script reads in .tif files of z-stack images, uses a weighted
% average of the intensity to find a surface, and uses this surface to
% calculate an indentation depth (and force).

clc, clear, close all

%%%% User Inputs
% Where are the images and where should output be saved?
directory = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\indentationDepth\imagesTIF\';
savedir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\indentationDepth\ExampleOutput\';

% Image resolution
umXY = 1.242; % um/pixel
umZ = 1; % um/pixel

% (Potentially) Resolution Dependent Image Analysis Parameters
gaussSmooth = 4; % Smooth noise in the image
cropSize = 140; % Size of region around indendation for further analysis
ringWidthFactor = 1.5; % Factor that remove reflections from the pipette from analysis
radiusRange = [40 70]; % In pixels, for finding the indentation

% Material Properties (for force calculations)
E = 0.2; % kPa;
vRange = [0.4 0.5]; % Range of Poisson values
R = 150/2; % um radius of pipette tip used for SLED

% Code Toggles
overwrite = 1; % Overwrite files if a .mat file for this image already exists? 0 = no, 1 = yes
slicePlot = 0; % Plot the found surface as a xz slice? 0 = no, 1 = yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Calculations:
%%%%% Set up data, import images, etc.
if ~exist(savedir,'file')
    mkdir(savedir)
end

list = dir([directory '*_channel1.tif']);

for imN = 1:length(list)
    
    imname = list(imN).name(1:end-4);
    disp(imname)
    
    if exist([savedir imname '_fitData.mat'],'file') && overwrite < 1
        disp('  Analysis already exists')
        continue
    end
    
    % Basic information about thte images
    N = imfinfo([directory list(imN).name]);
    sz(1) = N.Width;
    sz(2) = N.Height;
    N = length(N);

    % Load the z-stack
    im = NaN*ones(sz(1),sz(2),N);
    for jj = 1:N
        im(:,:,jj) = imread([directory list(imN).name],'tif',jj);
    end
    

    %%%%% Find the surface

    % Smooth image
    imSmooth = imgaussfilt(uint16(im),gaussSmooth);

    % Preallocate storage
    pos = NaN*ones(size(im,1),size(im,2));
    pos2 = pos;
    x = 1:size(im,3);
    for jj = 1:size(imSmooth,1)

        slice = double(squeeze(imSmooth(jj,:,:))');

        for kk = 1:size(slice,2)

            y = double(slice(:,kk));
            w = y/sum(y(:));
            pos(jj,kk) = sum(x.*w');

        end
        
        if slicePlot > 0

            figure(1)
            imagesc(squeeze(im(jj,:,:))')
            hold on
            plot(1:size(slice,2),pos(jj,:),'LineWidth',2,'Color','r')
            hold off
            set(gca,'DataAspectRatio',[1 umZ/umXY 1])
            title(jj)
            pause(0.01)    
            caxis([0 2^12])

            if mod(jj,100) == 0
               saveas(gcf,[savedir imname '_ExampleSliceSurface_y' num2str(jj) '.png'],'png')        
            end
        
        end

    end


    %%%% Scale to um and shift surface so bottom value is zero
    X = (1:size(im,1))*umXY;
    Y = (1:size(im,2))*umXY;
    posPlot = pos*umZ;
    posPlot = posPlot - min(posPlot(:)); % Bottom at zero

    %%%%%%%%%%% Figure of Heights
    figure(1)
    set(gcf,'Position',[100        306         560*1.5         420*1.5])

    surf(X,Y,posPlot)
    shading flat
    % set(gca,'DataAspectRatio',[1 1 umZ/umXY])
    set(gca,'DataAspectRatio',[1 1 1])
    axis([0 round(max(X(:))) 0 round(max(Y(:)))])
    set(gca,'FontSize',20)
    xlabel([char(181) 'm'],'FontSize',20)
    ylabel([char(181) 'm'],'FontSize',20)
    zlabel([char(181) 'm'],'FontSize',20)
    h = colorbar;
    set(get(h,'Label'),'String',['Height (' char(181) 'm)'])
    set(h,'FontSize',20)
    caxis([0 40])

    saveas(gcf,[savedir imname '_SledHeightMap.png'],'png')


    %%%% Create an ROI around the indentation    
    [centers, radii, metric] = imfindcircles(posPlot,radiusRange,'ObjectPolarity','Dark');
    if length(radii) > 1
        error('Two circles')
    elseif isempty(radii)
        error('No circles')
    end
    
    a = round(centers(2));
    b = round(centers(1));

    try
        indent = posPlot(a-cropSize:a+cropSize,:);
        indent = indent(:,b-cropSize:b+cropSize);
    catch
        warning('Small Box')
        a1 = max([1 a-cropSize]);
        a2 = min([a+cropSize size(posPlot,1)]);
        b1 = max([1 b-cropSize]);
        b2 = min([b+cropSize size(posPlot,2)]);
        indent = posPlot(a1:a2,:);
        indent = indent(:,b1:b2);
    end

    X2 = (1:size(indent,2))*umXY;
    Y2 = (1:size(indent,1))*umXY;

    %%%%%%%%%%% Figure of Indent as a Surface
    figure(2)
    surf(X2,Y2, indent)
    shading flat
    set(gca,'FontSize',20)
    xlabel([char(181) 'm'],'FontSize',20)
    ylabel([char(181) 'm'],'FontSize',20)
    zlabel([char(181) 'm'],'FontSize',20)

    saveas(gcf,[savedir imname '_ZoomIndent_NotToScale.png'],'png')
    
    
    %%%% Find the indendation location and depth
    % Make radial coordinates
    ring = ones(size(indent));
    [xr,yr] = meshgrid(1:size(ring,2),1:size(ring,1));
    
    % Find circles
    [centersI, radiiI, metricI] = imfindcircles(indent,radiusRange,'ObjectPolarity','Dark');
    if length(radiiI) > 1
        warning('Multiple Circles: Choosing largest radius (to keep the most data)')
        m = find(radiiI == max(radiiI));
        centersI = centersI(m,:);
        radiiI = radiiI(m);
    elseif isempty(radiiI)
        error('No Circles')
    end
    
    rr = sqrt((xr-centersI(1)).^2 + (yr-centersI(2)).^2);
    % Delete the ring around the indent
    ring(rr>radiiI & rr<ringWidthFactor*radiiI) = NaN;
    % Delete the reflection
    refl = ones(size(indent));
    refl( rr>ringWidthFactor*radiiI & xr>cropSize & yr>(centersI(2)-ringWidthFactor*radiiI) & yr<(centersI(2)+ringWidthFactor*radiiI) ) = NaN;
    
    % View the deleted version
    indentRing = indent.*ring.*refl;
    %%%%%%%%%%% Figure of Data with Reflection Removed
    figure(3)
    surf(indentRing)
    shading flat
    set(gca,'DataAspectRatio',[1 1 .5])
    view(2)
    set(gca,'FontSize',20)
    xlabel([char(181) 'm'],'FontSize',20)
    ylabel([char(181) 'm'],'FontSize',20)
    saveas(gcf,[savedir imname '_ZoomIndent_NotToScale_DataForFit.png'],'png')
   
    
    %%%%%% Calculate depth of indentation information    
    % baseline = area outside of the hole
    baseline = indentRing(rr>ringWidthFactor*radiiI);
    baseline(isnan(baseline)) = [];
    baseline = mean(baseline);
    % bottom = percentile to remove noise
    botIn = indentRing(rr<radiiI);
    botIn(isnan(botIn)) = []; % Shouldn't be any, but just in case of errors
    botIn = prctile(botIn,.1); % Min is probably noisy, use bottom hundredth percentile
    % depth = difference
    manDepth = baseline-botIn;
    disp(['     Depth (um): ' num2str(manDepth)])
    
    
    %%%%%%%%%%% Figure Showing depths
    figure(4)
    set(gcf,'Position',[100        306         560*1.5         420*1.5])
    % Downsample to make easier to view grid lines
    d = 2; % Downsampling factor
    X2d = downsample(X2,d);
    Y2d = downsample(Y2,d);
    indentRingd = downsample(indentRing,d);
    indentRingd = downsample(indentRingd',d)';
    % Surface of actual data
    s = surf(X2,Y2,indentRing);
    shading flat
    hold on
    % Mesh of downsampled data
    m = mesh(X2d,Y2d,indentRingd,'EdgeColor',0.2*[1 1 1],'FaceAlpha',0);
    hold off
    % Label axes
    xlabel([char(181) 'm'],'FontSize',20)
    ylabel([char(181) 'm'],'FontSize',20)
    zlabel([char(181) 'm'],'FontSize',20)
    grid on
    set(gca,'FontSize',20)
    view([-20.1000 9.5802])    
    %%% Add the baselines to the figure
    %%% (This might need tweaked to look 'pretty' for multiple replicates)
    hold on
    hbase = patch([0 400 400 0],[0 0 400 400],baseline*[1 1 1 1],0.8*[1 1 1],'EdgeColor',0.6*[1 1 1],'Facealpha',0.75);
    plot3(225*[1 1],10*[1 1],[botIn baseline],'-m','LineWidth',2)
    text(225+5,10,manDepth/2+botIn,[num2str(manDepth,'%0.1f') ' ' char(181) 'm'],...
        'Color','m','FontSize',20)
    hold off
    saveas(gcf,[savedir imname 'IndentLabeled_ShadingFacetedDownsampled.png'],'png')


    %%%% Herztian contact mechanics

    % F = 4/3 E/(1-v^2) d^3/2 R^1/2
    % v = Poisson's ratio
    % d = indentation depth
    % R = radius of pipette
    
    % Put in SI units
    dSI = manDepth/10^6;
    RSI = R/10^6; % m
    ESI = E*1000; % Pa = N/m^2
    
    Fpartial = 4/3*ESI*(dSI^1.5)*(RSI^.5); % in N divided by (1-v^2);

    Fmin = Fpartial/(1-vRange(1)^2);
    Fmax = Fpartial/(1-vRange(2)^2);

    disp(['     Force Range from Depth (nN): ' num2str(Fmin/10^-9) ' to ' num2str(Fmax/10^-9)])
    disp(' ')
    
    %%%% Save the data
    save([savedir imname '_fitData.mat'],...
        'umXY','umZ','gaussSmooth','cropSize','ringWidthFactor','radiusRange',...%input parameters
        'E','vRange','R',... %Material properites
        'posPlot','indent','X','Y','X2','Y2',...%height analysis
        'baseline','botIn','manDepth',...%manual depth analysis
        'd','Fmin','Fmax') % contact forces from manual depth 
    
end

%% Clean up
close all
disp('Batch Complete')