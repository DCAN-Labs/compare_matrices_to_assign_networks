# **Template Matching Documentation - supervised community detection for neural networks**

The human brain function is achieved through an emergent property of neurons to self-organize in a  repfuntion-feedback-reinforcement circuits.  Investigations into how these networks are represented on the brain, suggest that the typically-developing brain self-organizes into a general set of common networks. However, there is significant inter-subject variability in topography on the macroscopic scale (Cui et al. 2020; Glasser et al. 2016; Gordon, Laumann, Gilmore, et al. 2017; Gratton et al. 2018, 2020; Huth et al. 2016; Laumann et al. 2015; Rajkowska and Goldman-Rakic 1995; D. Wang et al. 2015; Brodmann 1909; von Economo and Koskinas 1925; Churchland and Sejnowski 1988).

This package  is designed to identfy neural networks using times series data. A series of analysis packages are used to quantify the networks topology afterwards.

The input data are assumed to be in processed with the freesurfer processing pipeline (theorhetically in BIDS format).  The code assumes that particpants have 91282 greyordinates  cortex + subcortical structures or are cortex only 549412).

The code works by taking taking a dense connectivity matrix and comparing the similarity to a series of network template from an independent data set.  If you don't have a dconn, you can build one with one of our other handy tools using the dense timeseries (dtseries.nii), see cifti_conn_matrix.m code for how to propoerly motion-censor the time series.


## **What is template matching?**
Template matching is a method of network mapping that leverages commonly observed networks that have been previously observed in a group average to accelerate the community detection process.  Individual networks are identified in a 3 step process:

1) **Run a community detection of your choosing on average connectivity matrix:**
In Hermosillo et al 2021, we used networks define from an average connectivit matrix created from 120 healthy adults. Infomap community detection was the conducted on the average matrix (See Gordon et al 2017 for details).

1) Generate a template of networks: Using the networks commmonly observed in a group 

3) Match connectivity of each voxel.


## Where is the code?

The template matching procesdure is part of a package called Compare matrices to Assign Networks.

The template matching package can be found at: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks

***

## Step-by-step tutorial



This package  is designed to identfy neural networks using times series data, specifically dense time series data in cifti format [(Glasser et al. 2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/). A series of analysis packages are used to quantify the networks topology afterwards.

The data are assumed to be in processed with the freesurfaer processing pipeline (theorhetically in BIDS format).  The code assumes that particpants have 91282 greyordinates  cortex + subcortical structures or are cortex only 549412).

The code works by taking taking a dense connectivity matrix and comparing the similarity to a series of network templates from an independet data set.  If you don't have a dconn, you can build one with one of our other handy tools using the dense timeseries (dtseries.nii), See [cifti_conn_matrix.m](https://github.com/DCAN-Labs/cifti-connectivity) code for how to properly motion-censor the time series.

%This code uses a connectivity matrices, a template connectivity matrix, and label file, to try to assing each greyordinate to a network for an individual.

## Step 0 : Gather the files your will need ahead of time.

### Required:
- **dconn** - the greyodrinate x time matrix where each cell contains the BOLD response.  It is highly recommended that this dconn is creaed using motion-regressed time series data in conjuction with motion scrubbing (If you're using the ABCD data set, you can use the .mat file within the derivatives that contains motion mask at various FD (framewise displacement thresholds). Click [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/) to read more about motion censoring and why it is important.
- **template file** - 
Each template file is different, 

There are several ways to do this:

1. If you have a correlation matrix for you subject and a template, this code will calculate the correlation between each a vector of correlations for greyordinate to every other greyodrinate.

Template matching.  This code will take a template of various networks and use caluculate an eta squared value against all other greyordinates.

## Step 1: Create your Template ... or download a premade one.

[![Placeholder Networks templates image](/assets/images/templates.png "Functional Network Templates"))


It is important to use an independent data set to identify the initial networks.  Using a starting set of networks, the code builds a seed-based correlation for all template subjects.

`makeCiftiTemplates_RH(dt_or_ptseries_conc_file,TR,all_motion_conc_file,project_dir,Zscore_regions,power_motion,remove_outliers,surface_only,use_only_subjects_that_pass_motion_criteria,combined_outliermask_provided)`

Alternatively, you can download our of our pre-made template from our repository: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/  

## Step 2: Make idividual-specific maps

The most current iteration of the code is *template_matching_RH.m*


`template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`
`


- subjectlist = subject (e.g. dconn.nii)
- data_type = "parcellated" or "dense" connectivity matrix
- template_path = path to .mat file that has th network templates.
- transform_data =  if you want to convert you data before comparing to your template, use can use 1 of 2 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearons' or 'Convert_to_Zscores' or use no tranformation
- output_cifti_name = output_cifti_name pretty clear
- cifti_output_folder = your project directory
- wb_command = path to run HCP workbench command.
- make_cifti_from_results = set to 1 if you want to save your results as a cifti. 0 will not save anything.
- allow_overlap = set to 1 if you're using overlapping networks in your cifti (Your input networks file will likely be a .dtseries.nii)
- overlap_method =  currently, the only supported method is "smooth_then_derivative"




#### Example call:

`template_matching_RH(/somepath/mydconn.nii, 'dense','somepathto/templatefile.mat',transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`



### Outputs are two files:
- eta squared value for each greyodrinate to specified network (.mat file).
- Second output file is the network assingment associated with maximum eta value.

## Step 3: 
