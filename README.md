# Template Matching Documentation 
### Supervised community detection for neural networks 


The human brain function is achieved through an emergent property of neurons to self-organize into neural networks.  Investigations into how these networks are represented on the brain, suggest that the typically-developing brain self-organizes into a general set of common networks. However, there is significant inter-subject variability in topography on the macroscopic scale (Cui et al. 2020; Glasser et al. 2016; Gordon, Laumann, Gilmore, et al. 2017; Gratton et al. 2018, 2020; Huth et al. 2016; Laumann et al. 2015; Rajkowska and Goldman-Rakic 1995; D. Wang et al. 2015; Brodmann 1909; von Economo and Koskinas 1925; Churchland and Sejnowski 1988).

This package  is designed to identify neural networks using times series data, specifically dense time series data in CIFTI format (dtseries.nii) [(Glasser et al. 2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/). A series of analysis packages are used to quantify the networks topology afterwards.


## **What is template matching?**
Template matching is a method of network mapping that leverages commonly observed networks that have been previously observed in a group average to accelerate the community detection process.  

Individual networks are identified in a 3-step process:

1) **Run a community detection of your choosing on average connectivity matrix:**
In Hermosillo et al. 2021, we used networks define from an average connectivity matrix created from 120 healthy adults. Infomap community detection was conducted on the average matrix (See Gordon et al. 2017 for details).

2) **Generate a template of networks:** 
Using the networks commonly observed in a group, we can create a template of networks in a specific file.

3) **Match connectivity of each voxel:** This will correlate each voxel, or seed, to its respective pair based on its function.

## Reasons for template matching

1. If you have a correlation matrix for your subject and a template file, this code will calculate the correlation between each a vector of correlations for greyordinate to every other greyodinate.

2. If you have a template of various networks, this code can use it to calculate an eta squared value against all other greyordinates.

## Where is the code?

The template matching procedure is part of a package called Compare Matrices to Assign Networks.

The template matching package can be found at: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks

***

## Step-by-step Tutorial: Background

The data are assumed to be processed with the FreeSurfer processing pipeline (theorhetically being in BIDS format). The code assumes that particpants have 91282 greyordinates  cortex + subcortical structures or are cortex only 549412.

The code works by taking taking a dense connectivity matrix (dconn.nii) and comparing its similarity to a series of network templates from an independent data set.  

If you don't have a dconn, you can build one with one of our other handy tools using the dense timeseries (dtseries.nii), See [cifti_conn_matrix.m](https://github.com/DCAN-Labs/cifti-connectivity) code for how to properly motion-censor the time series.

This code uses connectivity matrices, a template connectivity matrix, and a label file to try to assess each greyordinate to a network for an individual.

## Step 0: Gather the files you will need ahead of time

### Required:
- **dconn** - the greyodrinate x time matrix where each cell contains the BOLD response.  It is highly recommended that this dconn is creaed using motion-regressed time series data in conjuction with motion scrubbing (If you're using the ABCD data set, you can use the .mat file within the derivatives that contains motion mask at various FD (framewise displacement thresholds)). Click [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/) to read more about motion censoring and why it is important.
- **template file** - this will be a file of times series data (dtseries.nii) that is used to identify neural network connections.


If using SLURM, be sure to have ample time and space set out for the system to allocate the required resources, as the jobs can take up to 150 GB of RAM. 

## Step 1: Create your Template...or download a premade one

[![Placeholder Networks templates image](/assets/images/templates.png "Functional Network Templates")) 

This code is based on the 14 neural network-template defined by Gordon et al. 2017 seen above. (See [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub) for more information about each network.)


It is important to use an independent data set to identify the initial networks.  Using a starting set of networks, the code builds a seed-based correlation for all template subjects. This seed-based correlation shows how a seed voxel in one region of the brain is functionally related to another region of the brain based on the correlations between the time series of their activity.

`makeCiftiTemplates_RH(dt_or_ptseries_conc_file,TR,all_motion_conc_file,project_dir,Zscore_regions,power_motion,remove_outliers,surface_only,use_only_subjects_that_pass_motion_criteria,combined_outliermask_provided)`

Alternatively, you can download our of our pre-made template from our repository: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/  

## Step 2: Make idividual-specific maps

The most current iteration of the code is *template_matching_RH.m*


`template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`
`


- **dconn_filename =** path to input CIFTI file with network data 
- **data_type =** currently, the only supported data type is "dense"
- **template_path =** path to .mat file that has the network templates
- **transform_data =**  if you want to convert your data before comparing to your template, use can use 1 of 3 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearsons' or 'Convert_to_Zscores' or use no tranformation
    - "Covert_FisherZ_to_r" will convert your data through a Fisher-Z transformation (click [here](https://www.statisticshowto.com/fisher-z/) for more information)
    - "Convert_r_to_Pearsons" will convert your data into Pearson's coefficient (click [here](https://www.socscistatistics.com/tests/pearson/) for more information)
    - "Convert_to_Zscores" will convert your data into Z-scores (click [here](https://www.statisticshowto.com/probability-and-statistics/z-score/) for more information)
- **output_cifti_name =** name of the output CIFTI file
- **cifti_output_folder =** your project directory
- **wb_command =** path to run HCP workbench command
- **make_cifti_from_results =** "0" or "1"
    - "0" = does not save anything 
    - "1" = saves your results as a CIFTI file
-  **allow _overlap =** "0" or "1"
    - "0" = since your input networks file will likely be a .dtseries.nii, setting this to 0 will not use overlapping networks to analyze your data 
    - "1" = will use overlapping networks to analyze data
- **overlap_method =**  currently, the only supported method is "smooth_then_derivative"
- **surface_only =** "0" or "1" 
    - "0" = will not give output data on the surface
    - "1" = gives output data based on the surface 
- **already_surface_only =** "0" or "1"
    - "0" = use this if the input networks file does not have surface only data
    - "1" = use this if the input file is made up of surface data 



#### Example call:

`template_matching_RH(/somepath/mydconn.nii, 'dense','somepathto/templatefile.mat',transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`

[![Placeholder Example image])


### Outputs are two files:
1) A file with an eta squared value for each greyodrinate to specified network (.mat file).
2) A file with network assignments associated with maximum eta value (.dscalar file)

