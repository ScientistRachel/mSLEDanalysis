# calciumTimeLapse_noROI
This module of mSLEDanalysis measures fluorescent calcium dynamics in a cell monolayer in response to a uniform stimulus (e.g. H<sub>2</sub>O<sub>2</sub>).  These analyses are very similar to the analyses described in the calciumTimeLapse module, but do not use ROIs to subdivide the analysis.

## Dependencies
- The script `convertOIB_timeLapse.m` relies on OME's Bio-Format's toolbox for MATLAB to open the .oib microscope files.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).
- Some plots use the function `boundedline` to display error ranges: 
    > Kelly Kearney (2020). boundedline.m (https://www.github.com/kakearney/boundedline-pkg). Retrieved June 2, 2020.

## Running the Analyses
Required user inputs for each script are further described in each file.  The scripts should be run in this order:
1. `convertOIB_timeLapse.m` should be used to convert .oib microscope files into .tif files OR the user should make sure their images are saved as .tif files in the appropriate format.
2. `do_directory_setup.m` allows the user to specify groups of experiments (e.g. different cell lines).  The data structure set up in this script will be loaded in all following scripts to organize the data.
4. `find_dF_all_noROI.m` uses analyzes dF/F as a function of time.
5. Optional scripts to make additional plots or provide additional information are:
    * `plot_kymo_publicationFormat_noROI.m` creates kymographs with large fonts suitable for publication.
    * `plot_dF_all_noROI.m` combines multiple replicates into a single figure with error bars (requires `boundedline.m` and output from `find_dF_all_noROI.m`).
    * `report_values.m` reports summary statistics on dF/F to the command window (requires output from `find_dF_all_noROI.m`).

## Demo
Example images, as well as the expected output from the analyses, are available to [download at this link](https://drive.google.com/open?id=1HUxQJJmeYxVVFbY6DXMcV_50totPdLKF).
