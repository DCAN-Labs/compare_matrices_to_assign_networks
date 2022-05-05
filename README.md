# Template Matching Documentation 
### Supervised community detection for neural networks 

This package  is designed to identify neural networks using times series data, specifically dense time series data in CIFTI format (`.dtseries.nii`) [(Glasser et al. 2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/). A series of analysis packages are used to quantify the networks' topology afterwards. 

![Template Matching Method](/TM_README_Images/TM_RM_Image_1.jpg?raw=true) 


## What is template matching?
Template matching is a method of network mapping that leverages commonly observed networks that have been previously observed in a group average to accelerate the community detection process.

![14 Neural Networks](/TM_README_Images/TM_RM_Image_2.png?raw=true) 

This code uses a 14 neural network-base, first defined by [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub), and recreated with another group average. These 14 neural networks can be seen above. Read more about each network [here](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub).

Individual networks are identified in a 3-step process:

1) **Run a community detection of your choosing on average connectivity matrix:**
In Hermosillo et al. 2021, we used networks define from an average connectivity matrix created from 120 healthy adults. Infomap community detection was conducted on the average matrix (See [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub) for details).

2) **Generate a template of networks:** 
Using the networks commonly observed in a group, we can create a template of networks in a specific file.

3) **Match connectivity of each voxel:** This will correlate the connectivity of each voxel, or seed, to its respective pair, based on function.

## Reasons for template matching

1. If you have a correlation matrix for your subject and a template file, this code will calculate the correlation between each a vector of correlations for greyordinate to every other greyodinate.

2. If you have a template of networks, this code can use it to calculate the functional connectivity of each greyordinate in a subject with each template.
    - To explore probablistic network functions based on individual network maps for your data, check out [MIDB Atlas Maps](https://midbatlas.io).


## Background
The human brain function is achieved through an emergent property of neurons to self-organize into neural networks.  Investigations into how these networks are represented on the brain, suggest that the typically-developing brain self-organizes into a general set of common networks. However, there is significant inter-subject variability in topography on the macroscopic scale.

Refer to Cui et al. 2020; Glasser et al. 2016; Gordon, Laumann, Gilmore, et al. 2017; Gratton et al. 2018, 2020; Huth et al. 2016; Laumann et al. 2015; Rajkowska and Goldman-Rakic 1995; D. Wang et al. 2015; Brodmann 1909; von Economo and Koskinas 1925; Churchland and Sejnowski 1988 for more information. 


# The Code
## How do you install it?

To install the code, please download it from the following link: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/template_matching_RH.m

## How do you use it?
This code can only be run in the computer application **MATLAB**. You can download it here: https://www.mathworks.com/products/matlab.html?s_tid=mlh_so_learn. 

It can also be run from the command line using **MATLAB Runtime**. You can download that here: https://www.mathworks.com/products/compiler/matlab-runtime.html.

The MATLAB versions used here is **Version 9.6** (release name **R2019a**).

## What does this code do? 

This code uses connectivity matrices, a template connectivity matrix, and a label file to try to assess each greyordinate to a network for an individual connectivity matrix.


# Step-by-step Tutorial
## Background

A video tutorial for these steps can be found [here](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/TM_additional/TM_Video_Tutorial.mp4).

The data are assumed to be processed with the [FreeSurfer processing pipeline](https://github.com/DCAN-Labs/abcd-hcp-pipeline), preferably being in [BIDS format](https://bids-specification.readthedocs.io/). The code assumes that participants have 91282 greyordinates in the cortex and subcortical structures or have cortex-only data with 59412 greyordinates.

The code works by taking taking a dense connectivity matrix (`*.dconn.nii`) and comparing its similarity to a series of network templates (`*.dtseries.nii`) from an independent data set.  


## Step 1: Gather the files you will need ahead of time

### File 1. **Dense Connectivity Matrix** 

![Dense Connectivity Matrix Method and Example](/TM_README_Images/TM_RM_Image_3.png?raw=true) 

This will be a greyodrinate x time matrix in `.dconn.nii` format where each cell contains the BOLD response. 

These dense connectivity matrices are created by taking time series data of a single subject (A), analyzed to see which peaks (B) correspond to which neural network (C). The resulting data can be opened in Workbench View and altered, just as the dconn in (D) has been altered. 

- If you don't have a `*.dconn.nii`, you can build one with one of our other handy tools using the dense time series data (`*.dtseries.nii`) found [here](https://github.com/DCAN-Labs/abcd-hcp-pipeline) with the ABCD-HCP FreeSurfer processing pipeline. 

It is highly recommended that this file is created using motion-regressed time series data in conjuction with motion scrubbing. To properly motion-censor your data, refer to [cifti_conn_matrix.m](https://github.com/DCAN-Labs/cifti-connectivity). Click [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/) to read more about motion censoring and why it is important.

- If using the ABCD data set, you can use the `*.mat` file within the derivatives that contain motion mask at various framewise displacement (FD) thresholds. 


### File 2. **Template File**

![Template File Example](/TM_README_Images/TM_RM_Image_4.png?raw=true) 

This will be a file of times series data (`*.dtseries.nii`) that is used to identify neural network connections. It is important to use an independent data set to identify the initial networks. With a starting set of networks, the code builds a seed-based correlation for all template subjects. This seed-based correlation shows how a seed voxel in one region of the brain is functionally related to another region of the brain based on the correlations between the time series of their activity.

   - If you don't have a template file, proceed with one of two options:

     - **Option A: Create a template file**

        *Note:* Selection of subjects for template file should be considered carefully to ensure lack of bias. 

        Refer to [insert link to makeCiftiTemplates_README.md when finalized].

     - **Option B: Download a pre-made template file**

        Using a pre-made template file from the DCAN Lab's repository is also a viable option. It can be found here: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/support_files. 

*Important:* The time series data in the subject data set should correlate to the time series data in the template file for accurate matching. 

<br>
If using Slurm, be sure to have ample time and space set out for the system to allocate the required resources, as the jobs can take up to 150 GB of RAM. 

<br> 

## Step 2: Make individual-specific maps

The most current iteration of the code is **template_matching_RH.m** 


`template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`

A graphic with following inputs and outputs is added below. 

![Inputs and Outputs Graphic](/TM_README_Images/TM_RM_Image_5.jpg?raw=true) 


Below is a description of the required input parameters. 

- `dconn_filename` = path to input CIFTI file with network data 
- `data_type` = currently, the only supported data type is 'dense'
- `template_path` = path to .mat file that has the network templates
- `transform_data` =  if you want to convert the inputted data, use 1 of 3 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearsons' or 'Convert_to_Zscores'. Please note that the cortical and subcortical regions are transformed separately. You may also chose to use no tranformation, though the recommended transformation is 'Convert_to_Zscores'.
    - 'Covert_FisherZ_to_r' will convert your data through a Fisher-Z transformation (click [here](https://www.statisticshowto.com/fisher-z/) for more information)
    - 'Convert_r_to_Pearsons' will convert your data into Pearson's coefficient (click [here](https://www.socscistatistics.com/tests/pearson/) for more information)
    - 'Convert_to_Zscores' will convert your data into Z-scores (click [here](https://www.statisticshowto.com/probability-and-statistics/z-score/) for more information)
    - 'no_transformation' will not convert the input data. 
- `output_cifti_name` = name of the output CIFTI file
- `cifti_output_folder` = your project directory
- `wb_command` = [link](https://www.humanconnectome.org/software/connectome-workbench) to download Workbench Command to view outputs
- `make_cifti_from_results` = '0' or '1'
    - '0' = does not save anything 
    - '1' = saves your results as a CIFTI file
-  `allow _overlap` = '0' or '1'
    - '0' = since your input networks file will likely be a `*.dtseries.nii`, setting this to 0 will not use overlapping networks to analyze your data 
    - '1' = will use overlapping networks to analyze data
- `overlap_method` =  currently, the only supported method is 'smooth_then_derivative'
- `surface_only` = '0' or '1' 
    - '0' = will not give output data on the surface
    - '1' = gives output data based on the surface 
- `already_surface_only` = '0' or '1'
    - '0' = use this if the input networks file does not have surface only data
    - '1' = use this if the input file is made up of surface data 



#### Example call:

`template_matching_RH('/some/path/to/dconn_example.dconn.nii','dense','/some/path/to/template_example.mat','Convert_to_Zscores','example_run','/some/path/to/output/folder','/some/path/to/wb_command','1','1','smooth_then_derivative','0','0')`

*Note:* `/some/path/to/` refers to the computer path to the desired file. 

In this call:
- the subject data file is called **dconn_example** in `*.dconn.nii` format
- the data type is labeled **dense**
- the template file is called **template_example** in `*.mat` format
- the output data will be converted to Zscores
- the output data file will be named **example_run** in `*.dscalar` format
- the path to wb_command is written
- a CIFTI will be saved from the results as a `*.dscalar` file
- overlapping networks will be used to analyze data
- the overlap method will be **smooth_then_derivative**
- the code will generate outputs for all the data, not just surface data
- the data is known to not be surface only data 


# Outputs
There will be 1-5 outputs depending on the entered inputs.

![Intial MATLAB Outputs](/TM_README_Images/TM_RM_Image_6.png?raw=true) 

 1:. A `*.mat` file containing network names, new subject labels, and the functional connectivity of each greyordinate to the respective, provided template. 

 2: A `*.dscalar.nii` file with network assignments associated with maximum connectivity values. 
    
  -   3: This file will go through a cleaning script and generate a `*_recolored.dscalar.nii` file. 
    - *If requested,* this file will be saved locally.
    - *If only surface data is requested,* this file will contain the respective data.

 4: *If the overlapping networks are enabled*, an additional `*.dtseries.nii` file will be generated. 

 -  5: This file will also undergo a cleaning script, generating a `*_recolored.dtseries.nii` file.

# What can we do with these files? 

## Viewing the files

As some CIFTI files are outputted, you should view them with **Workbench View** (wb_view). These files will be the (1) `dscalar.nii`; (2) `_recolored.dscalar.nii`; (3) `dtseries.nii`; and (4) `_recolored.dtseries.nii`. You can download wb_view here: https://www.humanconnectome.org/software/connectome-workbench. It is part of a package by the Connectome Coordination Facility (CCF) called Connectome Workbench. 



An animated graphic showing the difference between the `dscalar.nii` and the `_recolored.dscalar.nii` files is shown below. 
![dscalars GIF](/TM_README_Images/TM_RM_Image_7.gif?raw=true) 

An animated graphic of resulting `dtseries.nii` data is shown below. This `dtseries.nii` file has been cleaned and is using the data from the `_recolored.dtseries.nii` file.
![dtseries GIF](/TM_README_Images/TM_RM_Image_8.gif?raw=true) 


## Further Analysis 

Using the resulting data of neural networks, you can calculate the perimeters of their borders using the following getborderperimeters.m code: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/getborderperimeters

You can also find the surface area of the resulting neural networks with the network_surface_area_from_network_file.m code following: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/network_surface_area