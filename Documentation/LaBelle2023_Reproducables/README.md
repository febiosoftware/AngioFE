# Directory Contents
This directory contains files to reproduce simulations from the pending manuscript "Spatial Configurations of 3D Extracellular Matrix Density and Anisotropy Simultaneously Guide Angiogenesis".
The subdirectories are detailed below:
* FEBio_AngioFE_Binaries
    * FEBio and AngioFE versions used for the publication.
* PseudoDeformationValidationMATLAB
    * MATLAB code used to validate pseudodeformation approach.
* Calibration_Simulations
    * Simulations of in vitro microvascular growth for calibrating simulation parameters.
* AnisotropyGradient_Simulations
    * Predictive simulations of microvascular growth in anisotropy gradients.
* TACS_Simulations
    * Predictive simulations of microvascular growth near tumor-associated collagen signatures (TACS)

# Reproducing AngioFE Simulations
Model inputs are provided as .feb files in the relevant directories.
#### FEBio and AngioFE distributions
Windows builds of FEBio and AngioFE are included in the FEBio_AngioFE_Binaries directory. Currently, AngioFE is compatible with FEBio 4.5 or FEBio 3. To run a model file, use the syntax:

`.\febio4.exe -noconfig -import .\AngioFE.dll -task=angio -i=[inputfile.feb]`

`.\febio3.exe -noconfig -import .\AngioFE.dll -task=angio -i=[inputfile.feb]`

#### Other distributions
Other distributions may be built from the appropriate FEBio and AngioFE commits. Input files are provided for AngioFE with FEBio 4.5. 
For the purpose of reproducibility, FEBio3 compatible versions are provided as well. The models for this publication were originally run with FEBio ver 3.7.c050b9169 and AngioFE ver 3.0.99e3b77 which can be built from:

https://github.com/febiosoftware/FEBio/commit/c050b9169

https://github.com/febiosoftware/AngioFE/commit/99e3b77


The syntax to run models from linux machines is:

`./febio4 -noconfig -import ./AngioFE.so -task=angio -i=[inputfile.feb]`

OR

`./febio3 -noconfig -import ./AngioFE.so -task=angio -i=[inputfile.feb]`

# Visualizing Results
#### Model outputs
Each model generates the following output files:
* FEBio .xplt file
    * The model results. This can be opened with FEBioStudio or later versions of PostView.
* AngioFE .ang2 vessel file
    *  A binary file that is used to visualize vessels within FEBioStudio.
*  FEBio .log file
    *  Log of FEBio outputs that are generated during model execution.
*  AngioFE _cells.txt file
    *  Text file detailing the position of tip cells at each model step.
*  AngioFE _log.csv file
    *  CSV file with microvessel network statistics for each time step.
    *  This file can be used to quickly determine the length and number of branches in the microvascular network at a given time.
*  AngioFE _time_stats.csv file
    *  CSV file with information on time spent during different functions.
*  AngioFE _vessels.csv file
    *  CSV containing information about the position and time that each segment originated at.
    *  This file can be used to determine microvessel network orientation distributions.

#### Visualizing results
1. Open the .xplt file in FEBioStudio.
2. Disable the mesh, disable the surface rendering, and change to orthographic projection (Windows key shortcuts m, w, and p).
3. Import the vessels via Post>Import Lines ...>Browse ... and select the corresponding .ang2 file.
4. Adjust line rendering. 
    1. Set the render mode to "smooth lines 3D"
    2. Set the line width to 8.