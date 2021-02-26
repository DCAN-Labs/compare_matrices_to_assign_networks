#! /bin/sh

# enter code directory
#cd "$( dirname "${BASH_SOURCE[0]}" )"

## Matlab command and usage
#mutualinfofromnetworks(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_vector,include_all_frames, smoothing_kernal,left_surface_file, right_surface_file, bit8, output_cifti_name,method, cifti_enhancement,other_half_networks_cii)
#e.g.mutualinfofromnetworks('/mnt/max/shared/projects/midnight_scan_club/data/MSC02_half1/MSC02_rfMRI_REST_FNL_preproc_Atlas.dtseries.nii','dtseries','/mnt/max/shared/projects/midnight_scan_club/data/MSC02_half1/power_2014_FD_only.mat', 0.2, 2.2,[1 2 3 4 5 10 15 20 25], 1 ,2.55,'/mnt/max/shared/projects/midnight_scan_club/data/MSC02_half1/MSC02.L.midthickness.32k_fs_LR.surf.gii','/mnt/max/shared/projects/midnight_scan_club/data/MSC02_half1/MSC02.R.midthickness.32k_fs_LR.surf.gii',0 , 'MSC02a_to_ADHD315_template_NE_MSC02b','template_matching',1,'MSC02_half2_matchedto_ADHD315_template_method_template_matching.dscalar.nii')
X="addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices')); mutualinfofromnetworks('${1}', '${2}', '${3}', '${4}', '${5}', '${6}' , '${7}', '${8}', '${9}', '${10}', '${11}', '${12}', '${13}', '${14}', '${15}', '${16}', '${17}', '${18}' )"

## This code then runs template matching or infomap on subjects and calculate the mututal information to a specifided dscalar.
##This code is designed to make correlation matrix from a subject's dtseries, and motion, and surface files (if smoothing is desired) See documentation for cifti_conn_matrix.
# Correlation matrices are automatically Z-scored.

#Hermosillo R. 1/25/2019
#--------------

#UPDATED 2/12/2019
#Updated to included option to run multiple repetitions 
#Update also one to sample the same dconn many times (in case you're interested in checking the stochatstic nature of infomap).

##How to run code
# Arguments are
# 1) dtseries = desnse time series
# 2) series = 'ptseries' or dtseries
# 3) motion_file = Power 2014 motion.mat file
# 4) FD_threshold = your frame-wise displament threshold. (e.g. the frames you want to keep from your time series.)
# 5) TR = repitition time of your data.
# 6) minutes_vector = a cell vector that contain each desired minute.  Leave empty to specifiy 'none minutes limit' (i.e. make an 'all frames' dconn)
# 7) include_all_frames, set to 1 to make a include minutes limit in your minutes vector
# 8) smoothing_kernal = specify smoothing kernal (note: if no smoothing file to be used, type 'none')
# 9) left_surface_file = full path to the left midthickness surface file.
# 10) right_surface_file full path to the right midthickness surface file.
# 11) bit8 = set to 1 if you want to make the outputs smaller in 8bit. Not recommended.  Set this to 0.
# 12) output_cifti_name = The name of your output cifti.  This option only works for template matching.  It does not work for infomap.
# 13) community_detection =  provide either 'template_matching' or 'Infomap'
# 14) method = Only used in template matching.  Specifiy 'Template matching'.
# 15) cifti_enhancement  = Not debugged.  Set to 0.
# 16) other_half_networks_cii,
# 17) num_interval_reps.  Number of repetitions of community detection. This is different that the specified number of reps that info map runs. Each repetition will make it's own dscalar.nii
# 18) indepen_time_series = if 1, each repititon will use a different dconn, if 0 it will use the same dconn. When minutes are specified these are randomly sampled from available frames below the FD threshold to make the dconn.  If you want to use the exact same dconn, this will make symlinks to that dconn. 

###########################################################################################
#matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2018a/matlab
matlab_exec=/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2016b/matlab

RandomHash=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16`
Tempmatlabcommand="matlab_command""$RandomHash"".m"

if [ -f "matlab_command""$RandomHash"".m" ]
then
	#echo "matlab_command.m found removing â¦"
	rm -fR "matlab_command""$RandomHash"".m"
fi

echo ${X} 
echo ${X} > "matlab_command""$RandomHash"".m"
cat "matlab_command""$RandomHash"".m"
${matlab_exec} -nodisplay -nosplash < "matlab_command""$RandomHash"".m"
rm -f "matlab_command""$RandomHash"".m"
