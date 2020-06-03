% Script to convert .oib files (microscope format) into .tif files for later
% analysis.  This script assumes that the .oib files have a single channel
% captured over multiple z-slices (one time point).
%
% Dependencies: This script requies the OME Bio-Formats MATLAB toolbox:
% https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html

clc, clear, close all

%%%%%%%%% User Inputs
directories = {'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\indentationDepth\imagesOIB\'}; % Folder of .oib files
savedirs = {'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\indentationDepth\imagesTIF\'}; % Location for .tif files
over_write = 0; % Overwrite files if a .tif file with that name already exists? 0 = no, 1 = yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Run Code

for mm = 1:length(directories)
    
    directory = directories{mm};
    list = dir([directory '*.oib']);
    
    savedir = savedirs{mm};
    if ~exist(savedir,'file')
        mkdir(savedir)
    end

    for kk = 1:length(list)
        
        if ~exist([savedir list(kk).name(1:end-4) '.tif',],'file') || over_write

            im = bfopen([directory list(kk).name]);
            im = im{1}; % Each row is an image; col 1 matrix & col 2 name

            for jj = 1:size(im,1)      
                if jj == 1
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '_channel1.tif'],'tif')
                else
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '_channel1.tif'],'tif','WriteMode','append')
                end
            end
    
        end

    end

end

clc
disp('Batch Complete')