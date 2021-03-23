# SOLID
Slicewise OutLier Detection for diffusion weighted MRI 

Tool has support for both Matlab and Python of which the Python version is currently most up to date. For details, please check the original paper [1] https://doi.org/10.1016/j.neuroimage.2018.07.003. Also, if you end up using this tool in your research, consider citing this article. The Matlab version can be used as an add-on to ExploreDTI (v.4.8.6) [2] with alternative subject motion and eddy current correction pipeline. The Python version can be used as an add-on to FSL's EDDY (6.0.1) [3] subject motion and eddy current correction pipeline. Both tools can be used as standalone to evaluate DWI quality.

----- Python -----

required arguments:  
  --in DWI         DWI file (*.nii.gz)  

optional arguments:  
  -h, --help       show this help message and exit  
  --bval BVAL      bval file (*.bval)  
  --bvec BVEC      bvec file (*.bvec)  
  --mask MASK      Mask file (*.nii.gz)  
  --eddyout EDDY   Path to precalculated FSL EDDY output  
  --out OUT        Output file  
  --metric METRIC  Used metric (var/mean/iod)  
  --PNG            Show modified Z-score as an image  
(  --GUI            Launch SOLID GUI) legacy  
  --thrU THRU      Set upper modified Z-score threshold (defaut 6.0)  
  --thrL THRL      Set lower modified Z-score threshold (defaut 3.5)  
  --smet SMET      Set scaling method for mapping modified Z-score to reliability weights: (linear) or sigmoid  
  --sfac SFAC      Set scaling factor for sigmoid scaling (defaut 0.2)  

A minimal example of quality control usage:  
> python SOLID.py --in data.nii.gz --bval data.bval --bvec data.bvec  
(Note: if file prefix is the same for all .nii, .bval, and .bvec the latter two can be ommited from input).

A minimal example how to implement in pipeline with FSL's EDDY:  
> eddy --imain=data.nii.gz ... --out=data_eddy_corrected  
> python SOLID.py --in data.nii.gz --eddyout data_eddy_corrected  
(Note: the argument for both EDDY --out and SOLID --eddyout must be the same string)  

----- Matlab -----
For a short YouTube tutorial about the Matlab tool, click the following figure: 

[![Alt text](https://img.youtube.com/vi/6R4tijOW4Ts/0.jpg)](https://www.youtube.com/watch?v=6R4tijOW4Ts "Click to show SOLID tutorial on youtube.com"). 

Matlab toolbox contains of three classes (tested on R2019a):  
@SOLID - all command line functionalities.  
@SOLID_GUI - improved graphical user interface. Requires @SOLID class.  
@SOLID_EDTI_plugin - a plugin for ExploreDTI v. 4.8.6 [2] Requires @SOLID and @SOLID_GUI.  

Extract these to your Matlab path or if you are using ExploreDTI to PathToExploreDTI/Source/. 

----- Notes and other tools -----

Modified Z-score maps (modZ4D) holds the outlier information in the original DWI space. This information can easily be summarized into 2D image where y-axis is the number of slices and x-axis is the number of volumes. However, after the motion correction transformation is applied to the 4D volume, these values much be considered as voxelwise information. Mapping of the modified Z-scores to weights can be achieved with various methods, of which linear and sigmoid are shown in the figure below:

![alt text](https://github.com/vilsaira/SOLID/blob/master/SOLID_modz2weight_mapping.png?raw=true)

Modified Z-score maps can be transformed using any registration algorithm using the transformation matrice obtained from the DWI motion and eddy correction steps. For example, Elastix (Transformix) [4] tool could be used like so  
 > transformix -in modZ4D.nii -tp TransformParameters.txt -out outputDirectory  

Only after this transformation to the corrected DWI space, the modified Z-scores should be mapped into reliability weights. By default, this mapping is done linearly between the set lower and upper thresholds. Typically, a modified Z-score of 3 is already 'outlierish' value and values above 6 tend to be extreme outliers but motion artefacts can produce modified Z-scores up to values of 20 therefore user can be quite liberal with the upper values. This likely depends on dataset. Python version also supports non-linear mapping and provides a sigmoid function for this step. By default sigmoid scaling factor is set to 0.2 which is somewhat between linear and full step function. By setting near zero scaling factor, user gets a step function. This can also be achieved by settings upper threshold close to the lower threshold with linear mapping.

How to use these reliability weights depends on the modelling software e.g. how diffusion tensor or some other model is fitted to the measurements. Reliability weights are intended to be used as a method to decrease the weight of unreliable data point in the modelling. The Matlab version of SOLID that works as ExploreDTI plugin is fully automated for WLLS and REKINDLE estimators. 

In microstructural modelling with MDT v0.14.5 [5] the reliability weights can be specified by the user as so:  
  > input_data = mdt.load_input_data(...volume_weights=weights)  

Many tools do not support robust modelling and these weights might not be easily added to pipelines usings such tools.  

----- References  -----
[1] Sairanen, V., A. Leemans, and C. M. W. Tax. "Fast and accurate Slicewise OutLIer Detection (SOLID) with informed model estimation for diffusion MRI data." NeuroImage 181 (2018): 331-346.  
[2] http://www.exploredti.com/  
[3] https://fsl.fmrib.ox.ac.uk/fsl/fslwiki  
[4] https://github.com/SuperElastix/elastix/wiki/FAQ  
[5] https://github.com/cbclab/MDT  
