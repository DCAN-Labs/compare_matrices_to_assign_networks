# Step-by-step Tutorial 
Written run through of template matching, the code of which can be found [here](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/template_matching_RH.m).


## 1. Acquire template file 
There are two ways to get a template file for this code. You can create a template file or download a premade one. First, a discussion on code-appropriate template files. 

### What is an appropriate template file?
Since this template matching tutorial uses the 14-neural network template as a basis for its matching, your template file should use the same divison for accurate matching (See [Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub) for details on the neural networks).

The code will use the template file, in accordance with a required subject data file (see step 2), to build a seed-based correlation matrix by showing functional relatedness between areas of the brain, or "voxels".

This template file should be a `.mat` file containing time series data. 

### Option A: Creating a template file

*Note:* Selection of subjects for template file should be considered carefully to ensure lack of bias. 

Using the code below, you can create a template file. 
`makeCiftiTemplates_RH(dt_or_ptseries_conc_file,TR,all_motion_conc_file,project_dir,Zscore_regions,power_motion,remove_outliers,surface_only,use_only_subjects_that_pass_motion_criteria,combined_outliermask_provided)` 

In this code, the input parameters refer to:
- `dt_or_ptseries_conc_file` = list of subjects' time series data; must be a `dtseries.nii` ending in `.conc`
- `TR` = frequency of every BOLD response in seconds
- `all_motion_conc_file` = list of subjects' motion data
    - for a high motion frame, it will be a `_motion.mat` file which is a standard output of DCAN processing pipeline
    - for a low motion frame, it will be a `_motion.txt` file  
- `project_dir` = path to output directory
- `Zscore_regions` = '0' or '1'
    - '0' = will be a Pearson correlation of the seed maps to the networks 
    - '1' = will be a Z-score transformation of the seed map to the networks
- `power_motion` = '0' or '1'
    - '0' = if using a `.txt file` for `all_motion_conc_file`
    - '1' = if using a `.mat file` for `all_motion_conc_file`
- `remove_outliers` = '0' or '1' 
    - '0' = will not remove outliers
    - '1' = will remove outliers 
- `surface_only` = '0' or '1' 
    - '0' = will not give output data on the surface data only
    - '1' = gives output data based on the surface data only
- `use_only_subjects_that_pass_motion_criteria` = '0' or '1' 
    - '0' = will use all provided subjects
    - '1' = will use the subjects that pass the defined motion criteria only  
- `combined_outliermask_provided` = '0' or '1' 
    - '0' = means a combined outlier mask is not provided
    - '1' = means a combined outlier mask is provided


### Option B: Downloading pre-made template file
Using a pre-made template file from the DCAN Lab's repository is also a viable option. It can be found here: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/tree/master/support_files. 
 

## 2. Add subject data set
To find the how areas in a certain subject's brain are related to each other in their function, you can use this code to match them according the neural network template discussed above.

### What is the subject data set?
This data set should be a `dconn.nii` file, a or dense connectivity matrix, consisting of greyordinates, a numerical map of the brain's structure, with respect to time. Each cell should contain the BOLD response. 

### How to build a dconn?

Refer to [`cifti_conn_wrapper.py`](https://github.com/DCAN-Labs/cifti-connectivity) code for how to properly motion-censor the time series. You can also use this code to create your `dconn.nii` file if you do not have one. 


*Important:* The time data in the subject data set should correlate to the template file time series data for accurate matching. 


## 3. Request system resources 
Running this code with the vastness of the subject data and the template file data will take a great deal of time and space. This applies to both pulling the resources and running the code if running on a cluster. 

For example: 
- Using SLURM to request resources can take up to 150 GB of RAM
- Using the code for all data can take up to 2 hours
    - On the other hand, using only surface data can take as little as 10 minutes 





## 4. Customize inputs 
For this code, there are many options of how you can analyze data including data transformations, making a CIFTI file from the results, using overlapping networks, and using only surface data to create the template. These inputs can be seen below. 

`template_matching_RH(dconn_filename, data_type,template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)`


Following is a description of the said inputs. 

- `dconn_filename` = path to input CIFTI file with network data 
- `data_type` = currently, the only supported data type is 'dense'
- `template_path` = path to .mat file that has the network templates
- `transform_data` =  if you want to convert your data before comparing to your template, use can use 1 of 3 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearsons' or 'Convert_to_Zscores' or use no tranformation
    - 'Covert_FisherZ_to_r' will convert your data through a Fisher-Z transformation (click [here](https://www.statisticshowto.com/fisher-z/) for more information)
    - 'Convert_r_to_Pearsons' will convert your data into Pearson's coefficient (click [here](https://www.socscistatistics.com/tests/pearson/) for more information)
    - 'Convert_to_Zscores' will convert your data into Z-scores (click [here](https://www.statisticshowto.com/probability-and-statistics/z-score/) for more information)
- `output_cifti_name` = name of the output CIFTI file
- `cifti_output_folder` = your project directory
- `wb_command` = path to run HCP workbench command
- `make_cifti_from_results` = '0' or '1'
    - '0' = does not save anything 
    - '1' = saves your results as a CIFTI file
-  `allow _overlap` = '0' or '1'
    - '0' = since your input networks file will likely be a .dtseries.nii, setting this to 0 will not use overlapping networks to analyze your data 
    - '1' = will use overlapping networks to analyze data
- `overlap_method` =  currently, the only supported method is 'smooth_then_derivative'
- `surface_only` = '0' or '1' 
    - '0' = will not give output data on the surface
    - '1' = gives output data based on the surface 
- `already_surface_only` = '0' or '1'
    - '0' = use this if the input networks file does not have surface only data
    - '1' = use this if the input file is made up of surface data 

## 5. Run the code

An example call is written below: 

`template_matching_RH('/some/path/to/dconn_example.dconn.nii','dense','/some/path/to/template_example.mat','Convert_to_Zscores','example_run','/some/path/to/output/folder','/some/path/to/wb_command','1','1','smooth_then_derivative','0','0')`

*Note:* `/some/path/to/` refers to the computer path to the desired file. 

In this call:
- the subject data file is called **dconn_example** in `.dconn.nii` format
- the data type is labeled **dense**
- the template file is called **template_example** in `.mat` format
- the output data will be converted to Zscores
- the output data file will be named **example_run** in `.dtseries.nii` format
- the path to wb_command is added
- a CIFTI will be saved from the results
- overlapping networks will be used to analyze data
- the overlap method will be **smooth_then_derivative**
- the code will generate outputs for all the data, not just surface data
- the data is known to not be surface only data 

Now, using MATLAB, run the code and wait patiently for your results.