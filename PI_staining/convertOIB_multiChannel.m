% Script to convert .oib files (microscope format) into .tif files for later
% analysis.  This script assumes that the .oib files have multiple channels
% but a single time point.
%
% Dependencies: This script requies the OME Bio-Formats MATLAB toolbox:
% https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html

clc, clear, close all

%%%%%%%%% User Inputs
directory = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\PI_staining\imagesOIB\'; % Folder of .oib files
savedir = 'D:\Code\_GitHubRepositories\SLEDanalysis\ExampleImages\PI_staining\imagesTIF\'; % Location for .tif files
over_write = 0; % Overwrite files if a .tif file with that name already exists? 0 = no, 1 = yes
channel_names = {'Calcium','PI','BrightField'}; % Channels in the order they are saved in the .oib file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Run Code
list = dir([directory '*.oib']);

for kk = 1:length(channel_names)
    if ~exist([savedir channel_names{kk} filesep],'file')
        mkdir(savedir,channel_names{kk})
    end
end

for kk = 1:length(list)
    
    im = bfopen([directory list(kk).name]);
    im = im{1}; % Each row is an image; col 1 matrix & col 2 name
    
    for jj = 1:size(im,1)
        
        imwrite(im{jj,1},[savedir channel_names{jj} filesep list(kk).name(1:end-4) '.tif'],'tif')
        
    end
    
end

clc
disp('Batch Complete')