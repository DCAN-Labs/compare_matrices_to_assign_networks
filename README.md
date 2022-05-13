# Template Matching Documentation 
### Supervised community detection for neural networks 

This package  is designed to identify neural networks using times series data, specifically dense time series data in CIFTI format (`.nii`) [(Glasser et al. 2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/). This method is diagrammed below. A series of analysis packages are used to quantify the networks' topology afterwards. 

![Template Matching Method](/TM_README_Images/template_matching_method.png?raw=true) 

(A): a template file's time series data has been overlaid with a dense connectivity matrix's neural network data, pointing out the DMN network as an example. (B): with the template matching code, individual network maps were created where the different networks are easily identifiable. 


## What is template matching?
Template matching is a method of network mapping that leverages commonly observed networks that have been previously observed in a group average to accelerate the community detection process.

![14 Neural Networks](/TM_README_Images/14_neural_networks.png?raw=true) 

This code uses a 14 neural network-base, first defined by [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub), and recreated with another group average seen above. Read more about each network [here](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub).

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

Refer to Cui et al. 2020; Glasser et al. 2016; Gordon, Laumann, Gilmore, et al. 2017; Gratton et al. 2018, 2020; Huth et al. 2016; Laumann et al. 2015; Rajkowska and Goldman-Rakic 1995; D. Wang et al. 2015; Brodmann 1909; von Economo and Koskinas 1925; and Churchland and Sejnowski 1988 for more information about the necessity of surface area data. 


# The Code
## How do you install it?

To install the code, please download it from the following link: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/template_matching_RH.m

## How do you use it?
This code can only be run in the computer application **MATLAB**. You can download it here: https://www.mathworks.com/products/matlab.html?s_tid=mlh_so_learn. 

It can also be run from the command line using **MATLAB Runtime**. You can download that here: https://www.mathworks.com/products/compiler/matlab-runtime.html.

The MATLAB versions used here is **Version 9.6** (release name **R2019a**).

## What does this code do? 

This code uses dense connectivity matrices, a template connectivity matrix, and a label file to try to assess each greyordinate to a network for its individual connectivity.


# Step-by-step Tutorial
## Background

A video tutorial for these steps can be found [here](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/TM_additional/TM_Video_Tutorial.mp4).

The data are assumed to be processed with the [FreeSurfer processing pipeline](https://github.com/DCAN-Labs/abcd-hcp-pipeline), preferably being in [BIDS format](https://bids-specification.readthedocs.io/). The code assumes that participants have 91282 greyordinates in the cortex and subcortical structures or have cortex-only data with 59412 greyordinates.

The code works by taking taking a dense connectivity matrix (`*.dconn.nii`) and comparing its similarity to a series of network templates (`*.dtseries.nii`) from an independent data set.  


## Step 1: Gather the files you will need ahead of time

### File 1. **Dense Connectivity Matrix** 

![Dense Connectivity Matrix Method and Example](/TM_README_Images/dconn_method_example.png?raw=true) 

This will be a greyodrinate x time matrix in `*.dconn.nii` format where each cell contains the [BOLD response](https://radiopaedia.org/articles/bold-imaging?lang=us). 

Dense connectivity matrices (dconns) are created by taking time series data of a single subject (A), analyzed to see which peaks (B) correspond to which neural network (C). The resulting data can be opened in Workbench View and altered, just as the resulting dense connectivity matrix in (D) has been altered. 

- If you don't have a `*.dconn.nii`, you can build one with one of our other handy tools using the dense time series data (`*.dtseries.nii`) found [here](https://github.com/DCAN-Labs/abcd-hcp-pipeline) with the ABCD-HCP FreeSurfer processing pipeline. 

It is highly recommended that this file is created using motion-regressed time series data in conjuction with motion scrubbing. To properly motion-censor your data, refer to [cifti_conn_matrix.m](https://github.com/DCAN-Labs/cifti-connectivity). Click [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/) to read more about motion censoring and why it is important.

- If using the ABCD data set, you can use the `*.mat` file within the derivatives that contain motion mask at various framewise displacement (FD) thresholds. 


### File 2. **Template File**

![Template File Example](/TM_README_Images/template_example.png?raw=true) 

This will be a file of times series data (`*.dtseries.nii`) that is used to identify neural network connections. In this template file, columns 4 and 6 correspond to neural networks that are determined to be inconsistent across subjects, and thus, omitted from the overall data.

It is important to use an independent data set to identify the initial networks. With a starting set of networks, the code builds a seed-based correlation for all template subjects. This seed-based correlation shows how a seed voxel in one region of the brain is functionally related to another region of the brain based on the correlations between the time series of their activity.

   - If you don't have a template file, proceed with one of two options:

     - **Option A: Create a template file**

        *Note:* Selection of subjects for template file should be considered carefully to ensure lack of bias. 

        Refer to makeCiftiTemplates_README.md for detailed instructions on how to create a template file.

     - **Option B: Download a pre-made template file**

        Using a pre-made template file from the DCAN Lab's repository is also a viable option. It can be found here: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/support_files. 

*Important:* The time series data in the subject data set should correlate to the time series data in the template file for accurate matching. 

<br>
If using Slurm, be sure to have ample time and space set out for the system to allocate the required resources, as the jobs can take up to 150 GB of RAM. 

<br> 

## Step 2: Make individual-specific maps

The most current iteration of the code is [**template_matching_RH.m**](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/template_matching_RH.m) 

As such, the current inputs will be: `template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`


Following is a description of the each of the required input parameters. 

- `dconn_filename` = path to `.dconn.nii` file with subject data
- `data_type` = currently, the only supported data type is **'dense'**
- `template_path` = path to `.mat` file that has the network templates
- `transform_data` =  if you want to convert the inputted data, use 1 of 3 transformations described below. 
    - **'Convert_to_Zscores'** will convert your data into Z-scores (click [here](https://www.statisticshowto.com/probability-and-statistics/z-score/) for more information)
    - **'Covert_FisherZ_to_r'** will convert your data through a Fisher-Z transformation (click [here](https://www.statisticshowto.com/fisher-z/) for more information)
    - **'Convert_r_to_Pearsons'** will convert your data into Pearson's coefficient (click [here](https://www.socscistatistics.com/tests/pearson/) for more information)
    - **'no_transformation'** will not convert the input data. 

         *Please note:* 
         
      - the cortical and subcortical regions are transformed separately
            
       - the recommended transformation is **'Convert_to_Zscores'**.
- `output_cifti_name` = name of the output CIFTI file(s)
- `cifti_output_folder` = the path to the project directory where you want to place these output files
- `wb_command` = [link](https://www.humanconnectome.org/software/connectome-workbench) to download Workbench Command to run code successfully
- `make_cifti_from_results` = '0' or '1'
    - '0' = does not save anything 
    - '1' = saves your results as a CIFTI file
-  `allow _overlap` = '0' or '1'
    - '0' = since your input networks file will likely be a `*.dtseries.nii`, setting this to 0 will not use overlapping networks to analyze your data 
    - '1' = will use overlapping networks to analyze data
- `overlap_method` =  currently, the only supported method is **'smooth_then_derivative'**
- `surface_only` = '0' or '1' 
    - '0' = will not give output data on the surface
    - '1' = gives output data based on the surface 
- `already_surface_only` = '0' or '1'
    - '0' = use this if the input networks file does not have surface only data
    - '1' = use this if the input file is made up of surface data 


#### Example call:

`template_matching_RH('/some/path/to/dconn_example.dconn.nii','dense','/some/path/to/template_example.mat','Convert_to_Zscores','testrun','/some/path/to/output/folder','/some/path/to/wb_command','1','1','smooth_then_derivative','0','0')`

*Note:* `/some/path/to/` refers to the computer path to the desired file. 

In this call:
- the subject data file is called **dconn_example** in `*.dconn.nii` format
- the data type is labeled **dense**
- the template file is called **template_example** in `*.mat` format
- the output data will be converted to **Zscores**
- the output data file will be named **testrun** in `*.dscalar` format
- the path to **wb_command** is written
- a CIFTI will be saved from the results as a `*.dscalar` file
- **overlapping networks** will be used to analyze data
- the overlap method will be **smooth_then_derivative**
- the code will generate outputs for **all data**, not just surface data
- the data is known to be **not surface only data** 

A variation of the example call's inputs are diagrammed below. 

![Inputs and Outputs Graphic](/TM_README_Images/inputs_outputs.png?raw=true) 

*Note:* For better image quality, please view the image at the following link: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/raw/master/TM_README_Images/inputs_outputs.png?raw=true.

As you can see, 5 different output files are generated. A detailed description of each is given in the following sections. 


# Outputs
There will be 1-5 outputs depending on the chosen inputs.


 1:. A `*.mat` file containing:
    
- (a) the functional connectivity scores of each greyordinate to the provided template
- (b) network names
- (c) new subject labels

 2: A `*.dscalar.nii` file with network assignments associated with maximum connectivity values. 
    
  -   3: This file will go through a cleaning script and generate a `*_recolored.dscalar.nii` file. 

      - *If requested,* this file will be saved locally.

        - *If only surface data is requested,* this file will contain the respective data.

 4: *If the overlapping networks are enabled*, an additional `*.dtseries.nii` file will be generated. 

 -  5: This file will also undergo a cleaning script, generating a `*_recolored.dtseries.nii` file.

# What can we do with these files? 

## Viewing the files

As the outputs include some CIFTI files (`.nii` files), you need to open them in **Workbench View** (wb_view).

You can download wb_view here: https://www.humanconnectome.org/software/connectome-workbench. It is part of a package by the Connectome Coordination Facility (CCF) called Connectome Workbench. 

Within wb_view, you will need blank brain image files on which you can overlay the outputs. A suggested directory of such templates is the [Midnight Scan Club (MSC)](https://github.com/MidnightScanClub/MSCcodebase). The GIFs below specifically use the "Conte69_atlas-v2.LR.32k_fs_LR.wb" files, found [here](https://github.com/MidnightScanClub/MSCcodebase/tree/master/Utilities/Conte69_atlas-v2.LR.32k_fs_LR.wb).

Once the brain templates have been loaded into wb_view, you can load in your CIFTI outputs. To view the network assignments, you can overlay the `*.dscalar.nii` file on the brain templates. 

Recommended settings for best viewing the neural networks is using the "power_surf" color palette, setting 17.85 as the positive maximum value  and 1.00 as the positive minimum value as shown below. 

![wb_view settings](/TM_README_Images/wb_view_settings.png?raw=true) 

As previously said, this `*.dscalar.nii` file undergoes a cleaning script, creating a `*_recolored.dscalar.nii` file. This cleaning script disregards connectivity values under a pre-set minimum value (shown in the `*.dscalar.nii` file as small, isolated "islands") as they are likely due to computer error and filters them to match whichever neighboring network the island shares most of its border with. The following animated graphic shows this difference between a `*.dscalar.nii` file and a `*_recolored.dscalar.nii` file.

![dscalars GIF](/TM_README_Images/difference_in_dscalars.gif?raw=true) 

Similarly, to view the `*.dtseries.nii` and the `*_recolored.dtseries.nii` files, you load in the brain templates and overlay the desired file(s). The display settings should be the same as when viewing the `*.dscalar.nii` files for consistency. As this file is with respect to time, meaning there are pictures at specific points in time, there are different results with different map settings. Below is an animated graphic showing how the results differ with map points 1-3. 

![dtseries different map points GIF](/TM_README_Images/dtseries_different_map_points.gif?raw=true) 

Like with the `*.dscalars.nii`, the `*.dtseries.nii` has undergone a cleaning script, creating the `*_recolored.dtseries.nii` file. In the graphic above, the `*_recolored.dtseries.nii` file has been overlaid. 

Following, all 16 map points of the `*_recolored.dtseries.nii` file are put together in an animated graphic. Map points 4 and 6 are empty as those networks are not used due to their observed inconsistency across subjects. Map point 14, corresponding to the Medial Temporal Lobe (MTL) neural network is empty only because it is not present in the pre-set left lateral view.  

![dtseries GIF](/TM_README_Images/dtseries_different_map_points.gif?raw=true) 



## Further Analysis 

Using the resulting data of neural networks, you can calculate the perimeters of their borders using the following getborderperimeters.m code: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/getborderperimeters

You can also find the surface area of the resulting neural networks with the network_surface_area_from_network_file.m code following: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/network_surface_area