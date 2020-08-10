# indentationDepth
This module of SLEDanalysis measures the indentation depth of the mSLED stimulus based on a fluorescent z-stack.  This indentation depth is further used to calculate an approximate applied force using Herztian contact mechanics.  The force calculation is approximate as it is based on the assumption that the cell monolayer's contribution to the deformation is negligible compared to the elastic modulus of the soft elastic dish.  This force calculation is based on the equation for a rigid sphere indenter in:
> C.T. McKee, J.A. Last, P. Russell, C.J. Murphy. _Indentation versus tensile measurements of Young's modulus for soft biological tissues_. Tissue Eng Part B Rev 2011, 17(3):155-164. [doi:10.1089/ten.teb.2010.0520]().

## Dependencies
The script `convertOIB_zstack.m` relies on OME's Bio-Format's toolbox for MATLAB to open the .oib microscope files.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).

## Running the Analyses
Required user inputs for each script are further described in each file.  The scripts should be run in this order:
1. `convertOIB_zstack.m` should be used to convert .oib microscope files into .tif files OR the user should make sure their images are saved as .tif files in the appropriate format.  This script expects the .oib file to contain multiple z-slices at a single time point for one fluorescent channels.
2. `measureIndentation.m` finds the weighted centroid of the fluorescence in z to create a surface.  The indentation is located by searching the image for circles, and the depth is calculated as the distance from the bottom of the indentation to the average baseline height.

## Demo
Example images, as well as the expected output from the analyses, are available to [download at this link](https://drive.google.com/open?id=1NBGhNU0JM8EceidWO0Z7DkPLCvT_APmY).
