# SOLID
Slicewise OutLier Detection for diffusion weighted MRI  

For details, please check the original paper [1] https://doi.org/10.1016/j.neuroimage.2018.07.003. For a short YouTube tutorial, click the following figure: 

![Alt text](https://img.youtube.com/vi/6R4tijOW4Ts/0.jpg)](https://www.youtube.com/watch?v=6R4tijOW4Ts "Click to show SOLID tutorial on youtube.com"). 

Matlab toolbox contains of three classes (tested on R2019a):
@SOLID - all command line functionalities. 
@SOLID_GUI - improved graphical user interface. Requires @SOLID class.
@SOLID_EDTI_plugin - a plugin for ExploreDTI v. 4.8.6 [2] Requires @SOLID and @SOLID_GUI.

Extract these to your Matlab path or if you are using ExploreDTI to PathToExploreDTI/Source/. 

Python script can be used for bulk investigation of outliers if Matlab is not available. Help files are included within the source e.g. in Matlab type ">help SOLID" and in python "python SOLID.py -h".

1.0) In general, to calculate modified Z-scores i.e. find outliers
--> Can use Matlab script, Matlab GUI, or Python script for this

2.0) Transformation of modified Z-score maps after motion/distortion correction

2.1) Commandline example for Elastix (Transformix) [3]
 > transformix -in modZ.nii -tp TransformParameters.txt -out outputDirectory

Prerequisites this are to install Elastix tools to OS path and use them or any other compatible tool to perform image registration on DWIs to obtain TransformParameters file. These same transformation are simply applied to the modZ nifti.

2.2) With FSL [4] in progress...

3.0) Apply thresholding according to eq. 4 in [1] on the transformed modified Z-scores (i.e. output of step 2).

4.0) Use SOLID weights in the model estimation. This depends highly on the software you are using. For example, MDT v0.14.5 and later has the required functionality to add user defined weights [5]:

  > input_data = mdt.load_input_data(...volume_weights=weights)
                                
References:
[1] Sairanen, V., A. Leemans, and C. M. W. Tax. "Fast and accurate Slicewise OutLIer Detection (SOLID) with informed model estimation for diffusion MRI data." NeuroImage 181 (2018): 331-346.
[2] https://github.com/SuperElastix/elastix/wiki/FAQ
[3] https://fsl.fmrib.ox.ac.uk/fsl/fslwiki
[4] https://github.com/cbclab/MDT
[5] http://www.exploredti.com/
