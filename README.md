# Analyzing *M*echanical *S*timulation on *L*ow *E*lastic Modulus *D*ishes (mSLED)

This code was developed to measure responses of a cell monolayer to mSLED.  The code is separated into modules which can be run separately as desired.  The scripts save .mat files and image files of figures as output.  Those .mat files can be used to further compare the output of the different modules. For example, the dF/F measurements from _calciumTimeLapse_ can be compared to the cell death measurements from _PI_staining_.  These analyses were used to investigate mechanotransduction in breast cancer cells by Pratt, et al., which should be cited in further use of this code:
> S.J. Pratt , R. M. Lee , K. T. Chang , E. O. Hernández-Ochoa , D. A. Annis , E. C. Ory , K. N. Thompson , P. C. Bailey , T. J. Mathias , J. A. Ju , M. I. Vitolo , M. F. Schneider , J. P. Stains , C. W. Ward , and S. S. Martin. _Mechanoactivation of NOX2-generated ROS elicits persistent TRPM8 Ca<sup>2+</sup> signals that are inhibited by oncogenic KRas_. 2020.

## Dependencies
All scripts were developed and tested in MATLAB 2019a with the Image Processing Toolbox.  Some scripts further require additional dependencies:
- The scripts `convertOIB_*.m` rely on OME's Bio-Format's toolbox for MATLAB to open the .oib microscope files.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).
- Some plots use the function `boundedline` to display error ranges: 
    > Kelly Kearney (2020). boundedline.m (https://www.github.com/kakearney/boundedline-pkg). Retrieved June 2, 2020.

## Modules
Instructions for running each module and additional details about the code are available in each subfolder's README.
- _PI_staining_ measures changes in propidium iodide (PI) staining across two images (e.g. before and after scratching).
- _calciumTimeLapse_ measures dF/F as a function of time and as a function of distance from the scratch.
- _calciumTimeLapse_noROI_ measures dF/F as a function of time for uniform stimulations, such as H<sub>2</sub>O<sub>2</sub>.
- _indentationDepth_ measures the depth of the scratch and resulting applied force on a soft elastic dish.

## Demo
Example images, as well as the expected output from the analyses, are available to [download at this link](https://drive.google.com/drive/folders/1XwqTBv__Ahuj9axImJj1oEZcCsMgein-?usp=sharing).  The user should change the image directory at the top of the relevant analysis script to match the downloaded example images and then run the script without changing other parameters to create the same example output.
