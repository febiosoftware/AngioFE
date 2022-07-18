AngioFE code changes
======================

1. SetParamAttribute is no longer part of febio4. Need alternative mechanism for processing these parameters. (Currently, commented out)

3. Made changes to base classes (mostly replaced FEMaterial with FEMaterialProperty). 

4. Made a change to how mat_axis is processed. 


TODO:

1. Global constants are not discoverable by febio studio 2. Implement alternative approach of defining global variables (e.g. make them task parameters)

5. Files need to be converted to febio3 or febio4 format. 
