# PI_staining
This module of mSLEDanalysis measures changes in propidium iodide (PI) staining across two images (e.g. images before and after scratching).  It is assumed that there is no translation of the sample between images.  If the sample was moved between imaging, the images would need to be registered before further analysis.  Images must be captured with the same acquisition parameters for this analysis to be meaningful.

## Dependencies
The script `convertOIB_multiChannel.m` relies on OME's Bio-Format's toolbox for MATLAB to open the .oib microscope files.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).

## Running the Analyses
Required user inputs for each script are further described in each file.  The scripts should be run in this order:
1. `convertOIB_multiChannel.m` should be used to convert .oib microscope files into separate .tif files OR the user should make sure their images are saved as .tif files in the appropriate format.  This script expects the .oib file to contain a single z-slice at a single time point for multiple imaging channels.
2. `PI_staining_analysis` compares paired images (e.g. "Prescratch" and "Postscratch") and calculates the difference in PI staining intensities.

## Demo
Example images, as well as the expected output from the analyses, are available to [download at this link](https://drive.google.com/open?id=18n_1H27cwbsWqO9va_ekhjuvqBxfZHCE).
