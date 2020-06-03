% Script to convert .oib files (microscope format) into .tif files for later
% analysis.  This script assumes that the .oib files have two channels
% captured over time.
%
% Dependencies: This script requies the OME Bio-Formats MATLAB toolbox:
% https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html

clc, clear, close all

%%%%%%%%% User Inputs
directories = {'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\imagesOIB\'}; % Folder of .oib files
savedirs = {'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\calciumTimeLapse_noROI\imagesTIF\'}; % Location for .tif files
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

            im = bfopen([directory list(kk).name]); % This line requires the bioformats 
            im = im{1}; % Each row is an image; col 1 matrix & col 2 name

            % Other format has fluorescent images first then brightfield, so
            % format the same way for consistency
            for jj = 1:2:size(im,1)   % Calcium         
                if jj == 1
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '.tif'],'tif')
                else
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '.tif'],'tif','WriteMode','append')
                end
            end
            for jj = 2:2:size(im,1)  % Bright field
                if jj == 1
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '.tif'],'tif')
                else
                    imwrite(im{jj,1},[savedir list(kk).name(1:end-4) '.tif'],'tif','WriteMode','append')
                end
            end
        
        end

    end

end

clc
disp('Batch Complete')