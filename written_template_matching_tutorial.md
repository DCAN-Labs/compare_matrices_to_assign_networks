# Step-by-step Tutorial 
Written run through of template matching, the code of which can be found [here](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/template_matching_RH.m). 

## 1. Acquire template file 
There are two ways to get a template file for this code. You can create a template file or download a premade one. First, a discussion on code-appropriate template files. 

### What is an appropriate template file?
Since this template matching tutorial uses the [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub) 14-neural network template as a basis for its matching, your template file should use the same divison for accurate matching. 

The code will use the template file, in accordance with a required subject data file discussed below (See Step 2), to build a seed-based correlation matrix by showing functional relatedness between areas of the brain, or "voxels".

This template file should be a *.mat file containing time series data. 

### Option A: Creating a template file
Using the code below, you can create a template file. 
`makeCiftiTemplates_RH(dt_or_ptseries_conc_file,TR,all_motion_conc_file,project_dir,Zscore_regions,power_motion,remove_outliers,surface_only,use_only_subjects_that_pass_motion_criteria,combined_outliermask_provided)` 

### Option B: Downloading pre-made template file
Using a pre-made template file from the FairLab's repository is also a viable option, found here: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/. 
 

## 2. Add subject data set
To find the how areas in a certain subject's brain are related to each other in their function, you can use this code to match them according the neural network template discussed above.

### What is the subject data set?
This data set should be a dconn.nii file, a or dense connectivity matrix, consisting of greyordinates, a numerical map of the brain's structure, with respect to time. Each cell should contain the BOLD response. 

*Important:* The time data in the subject data set should correlate to the template file time series data for accurate matching. 


## 3. Request system resources 
Running this code with the vastness of the subject data and the template file data and the usage of a workbench command will take a great deal of time and space. This applied to both pulling the resources and running the code. For example, using SLURM to request resources can take up to 150 GB of RAM. 

## 4. Customize inputs 
For this code, there are many options of how you can analyze data including data transformations, making a CIFTI file from the results, using overlapping networks, and using only surface data to create the template. These inputs can be seen below. 

`template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`

For a more detailed explanation of each input, refer to the [README.md](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/README.md) file for this code. 

## 5. Run the code
Using MATLAB, run the code and wait patiently for your results. 

An example call is written below: 

`template_matching_RH(/somepath/mydconn.nii, 'dense','somepathto/templatefile.mat',transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`
