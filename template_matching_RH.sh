#! /bin/sh

# enter code directory
#cd "$( dirname "${BASH_SOURCE[0]}" )"

## Matlab command and usage

#template_matching_RH(subjectlist, data_type, template_path,transform_data,output_cifti_name, wb_command)
#e.g. template_matching_RH('/home/exacloud/lustre1/fnl_lab/data/HCP/processed/midnight_scan_club/MSC01/20180509-date/HCP_release_20170910_v1.3/MSC01/MNINonLinear/Results/MSC01_rfMRI_REST_FNL_preproc_Atlas_SMOOTHED_2.55.dtseries.nii_all_frames_at_FD_0.2.dconn.nii','dense','/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_fromsurfonly.mat','Convert_FisherZ_to_r','MSC01_to_ADHD315','/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-9253ac2/bin_rh_linux64/wb_command')

#twins_mapping_wrapper(dt_or_ptseries_conc_file(1),motion_file(2),left_surface_file(3), right_surface_file(4), output_file_name(5), cifti_output_folder(6),TR(7),minutes_limit(8),FD_threshold(9),transform_data(10),template_path (11),surface_only(12),already_surface_only(13), use_all_ABCD_tasks (14), run_infomap_too (15) ,output_directory (16), dtseries_conc (17), use_continous_minutes (18), memory_limit_value (19),clean_up_intermediate_files (20), wb_command (21), additional_mask (22)), remove_outliers(23)

X="addpath('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks'); template_matching_RH('${1}', '${2}', '${3}','${4}', '${5}', '${6}', '${7}', '${8}', '${9}', '${10}', '${11}', '${12}')"
#function twins_mapping_wrapper(dt_or_ptseries_conc_file,motion_file,left_surface_file, right_surface_file, output_file_name, cifti_output_folder)
#function twins_mapping_wrapper(dt_or_ptseries_conc_file,motion_file,left_surface_file, right_surface_file, output_file_name, cifti_output_folder,TR,minutes_limit,FD_threshold)
#Hermosillo R. 4/19/2019
#this code runs template matching starting from a dtseries.  Several parameters are hardcoded into the corresponding matlab code.


#$1= path to dconn.nii (or pconn) file

###########################################################################################
#matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2018a/matlab
#matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2016b/matlab
#matlab_exec=/panfs/roc/msisoft/matlab/R2019a/bin/matlab
matlab_exec=/common/software/install/migrated/matlab/R2019a/bin/matlab
RandomHash=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16`
Tempmatlabcommand="matlab_command""$RandomHash"".m"

if [ -f "matlab_command""$RandomHash"".m" ]
then
	#echo "matlab_command.m found removing â¦"
	rm -fR "matlab_command""$RandomHash"".m"
fi

#echo ${X} 
echo ${X} > "matlab_command""$RandomHash"".m"
cat "matlab_command""$RandomHash"".m"
${matlab_exec} -nodisplay -nosplash < "matlab_command""$RandomHash"".m"
rm -f "matlab_command""$RandomHash"".m"
