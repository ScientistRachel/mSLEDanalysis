# calciumTimeLapse
This module of SLEDanalysis measures fluorescent calcium dynamics in a cell monolayer in response to scratching.  These analyses were originally described in a publication by Pratt, et al:
> S. J. Pratt, E. O. Hernández-Ochoa, R. M. Lee, E. C. Ory, J. S. Lyons, H. C. Joca, A. Johnson, K. Thompson, P. Bailey, C. J. Lee, T. Mathias, M. I. Vitolo, M. Trudeau, J. P. Stains, C. W. Ward, M. F. Schneider, and S. S. Martin. _Real-time scratch assay reveals mechanisms of early calcium signaling in breast cancer cells in response to wounding_. Oncotarget, 9:25008-25024, 2018. [doi:10.18632/oncotarget.25186](http://dx.doi.org/10.18632/oncotarget.25186)

## Dependencies
- The script `convertOIB_timeLapse.m` relies on OME's Bio-Format's toolbox for MATLAB to open the .oib microscope files.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).
- Some plots use the function `boundedline` to display error ranges: 
    > Kelly Kearney (2020). boundedline.m (https://www.github.com/kakearney/boundedline-pkg). Retrieved June 2, 2020.

## Running the Analyses
Required user inputs for each script are further described in each file.  The scripts should be run in this order:
1. `convertOIB_timeLapse.m` should be used to convert .oib microscope files into .tif files OR the user should make sure their images are saved as .tif files in the appropriate format.
2. `do_directory_setup.m` allows the user to specify groups of experiments (e.g. different cell lines).  The data structure set up in this script will be loaded in all following scripts to organize the data.
3. `find_tip_all.m` uses the brightfield microscopy images to find the location of the scratching pipette.  This location will be used to find ROIs around the SLED region and analyze metrics as a function of distance from the scratch.
4. The dF/F analyses can be run separately and in any order:
    * `find_dF_all_v2.m` uses ROIs to analyze dF/F as a function of time.
    * `find_dF_dist_all_v2.m` analyzes dF/F as a function of distance from the scratch.
5. Optional scripts to make additional plots or provide additional information are:
    * `plot_kymo_publicationFormat.m` creates kymographs with large fonts suitable for publication (requires output from `find_tip_all.m`).
    * `plot_dF_all_v2.m` combines multiple replicates into a single figure with error bars (requires `boundedline.m` and output from `find_dF_all_v2.m`).
    * `plot_dF_dist_all_v3.m` combines multiple replicates into a single figure with error bars (requires `boundedline.m` and output from `find_dF_dist_all_v2.m`).
    * `report_values.m` reports summary statistics on dF/F to the command window (requires output from `find_dF_all_v2.m` and `find_dF_dist_all_v2.m`).

## Demo
Example images, as well as the expected output from the analyses, are available to [download at this link](https://drive.google.com/open?id=1vJPFLyV_lPtSy1lvhVsndyrbKgiUpFOa).
