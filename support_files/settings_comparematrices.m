function settings=settings_comparematrices()
%% Define the paths for support functions
str = computer;
code_dir = pwd;

try
    disp(['code dir is : ' code_dir]);
    server_name=code_dir(1:10);
catch
    disp(['code dir is : ' code_dir]);
    try
        server_name=code_dir(1:5);
    catch
     disp(['code dir is : ' code_dir]);       
        server_name=code_dir(1:4); % assume the server is /tmp
    end
end


switch server_name
    
    case '/code' %singularity container
        path{1}='/code/support_files/gifti-1.6'; % added during compiling
        path{2}='/code/support_files/Matlab_CIFTI'; % added during compiling
        path{3}='/code/support_files/Matlab_effect_size_toolbox/'; % added during compiling
        path{4}='/code/support_files/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii'; %unused.
        %templates:
        path{5}='/code/support_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/code/support_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/code/support_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii';
        path{8}='/code/support_files/91282_Greyordinates.dscalar.nii';
        path{9}='/code/support_files/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii'; %unused
        path{10}='/code/support_files/91282_Greyordinates.dtseries.nii';
        path{11}='/code/support_files/91282_Greyordinates_surf_only.dtseries.nii'; 
        path{12}='/code/support_files/91282_Greyordinates_surf_only.dscalar.nii';              
        path{13}='/code/support_files/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene';      
        path{14}='/code/support_files/MSC01_template_scene_subcort_label_MSI.scene';      
        path{15}='/code/support_files/make_dscalar_pics_v9.3.sh';
        path{16}='/code/support_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
        
        path_wb_c='/opt/workbench/wb_command'; % workbench command path
        path_template_nets='/code/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat'; %unused - supplied by user in template_matching_RH.

    case '/mnt/max/s' %rushmore
        path{1}='/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/';
        path{2}='/mnt/max/shared/code/external/utilities/Matlab_CIFTI';
        path{3}='/mnt/max/shared/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/mnt/max/shared/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/10050-1/20100430-SIEMENS_TrioTim-Nagel_K_Study/HCP_release_20161027/10050-1/MNINonLinear/Results/10050-1_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii_to_Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii.pscalar.nii';
        path{8}='/mnt/max/shared/code/internal/pipelines/HCP_release_20170910_v1.4/global/templates/91282_Greyordinates/91282_Greyordinates.dscalar.nii'; 
        path{9}='/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
        path{10}='/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii';              
        path{11}='/mnt/max/shared/code/external/utilities/gifti-1.6';
        path{12}='/mnt/max/shared/code/internal/utilities/Zscore_dconn/';
        path{13}='/mnt/max/shared/code/internal/analyses/compare_matrices/support_files/91282_Greyordinates.dtseries.nii';
        %path{12}='/mnt/max/shared/code/internal/analyses/compare_martices/support_files/seedmaps_ADHD_Control_20_subs_dtseries_all_networks_fromsurfonly.mat';
        
        path_template_nets='/mnt/max/shared/code/internal/analyses/compare_matrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_fromsurfonly.mat';
        path_wb_c='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command'; % workbench command path
  
    case '/home/exac' %lustre filesystem
        path{1}='/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities';
        path{2}='/home/exacloud/lustre1/fnl_lab/code/external/utilities/Matlab_CIFTI';
        path{3}='/home/exacloud/lustre1/fnl_lab/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/home/exacloud/lustre1/fnl_lab/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/10050-1/20100430-SIEMENS_TrioTim-Nagel_K_Study/HCP_release_20161027/10050-1/MNINonLinear/Results/10050-1_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii_to_Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii.pscalar.nii';
        path{8}='/home/exacloud/lustre1/fnl_lab/code/internal/pipelines/HCP_release_20170910_v1.4/global/templates/91282_Greyordinates/91282_Greyordinates.dscalar.nii'; 
        path{9}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
        path{10}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii';
        path{11}='/home/exacloud/lustre1/fnl_lab/code/external/utilities/gifti-1.6';
        path{12}='/home/exacloud/lustre1/fnl_lab/code/internal/utilities/Zscore_dconn/';
        path{13}='/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/91282_Greyordinates.dtseries.nii';
        path_template_nets='/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';       
        path_wb_c='/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-1.2.3-HCP/bin_rh_linux64/wb_command';
          
    case '/home/fair' %mesabi
        path{1}='/home/faird/shared/code/external/utilities/gifti-1.6';
        path{2}='/home/faird/shared/code/internal/utilities/Matlab_CIFTI';
        path{3}='/home/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/home/faird/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/home/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii';
        path{8}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii'; 
        path{9}='/home/faird/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';        
        path{10}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dtseries.nii';
        path{11}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dtseries.nii';
        path{12}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii';      
        path{13}='/home/faird/shared/code/internal/utilities/figure_maker/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene';      
        path{14}='/home/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_scalar_MSI.scene'; % avoid using the following, since it is intended for a label file instead of dscalar file. '/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_label_MSI.scene';      
        path{15}='/home/faird/shared/code/internal/utilities/figure_maker/make_dscalar_pics_v9.4.sh';
        path{16}='/home/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
        
        path_wb_c='/home/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_template_nets='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';
        
    case '/panfs/roc' %mesabi
        path{1}='/panfs/jay/groups/6/faird/shared/code/external/utilities/gifti-1.6';
        path{2}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/Matlab_CIFTI';
        path{3}='/panfs/jay/groups/6/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii';
        path{8}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii';
        path{9}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
        path{10}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dtseries.nii';
        path{11}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dtseries.nii'; 
        path{12}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii';              
        path{13}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene';      
        path{14}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_scalar_MSI.scene'; % avoid using the following, since it is intended for a label file instead of dscalar file. '/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_label_MSI.scene';      
        path{15}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/make_dscalar_pics_v9.4.sh';
        path{16}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
        
        %path_wb_c='/panfs/roc/groups/4/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_wb_c='/panfs/jay/groups/2/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_template_nets='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';

    case '/panfs/jay' %mesabi
        path{1}='/panfs/jay/groups/6/faird/shared/code/external/utilities/gifti-1.6';
        path{2}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/Matlab_CIFTI';
        path{3}='/panfs/jay/groups/6/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii';
        path{8}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii';
        path{9}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
        path{10}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dtseries.nii';
        path{11}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dtseries.nii'; 
        path{12}='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii';              
        path{13}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene';      
        path{14}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_scalar_MSI.scene'; % avoid using the following, since it is intended for a label file instead of dscalar file. '/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_label_MSI.scene';      
        path{15}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/make_dscalar_pics_v9.4.sh';
        path{16}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
        
        %path_wb_c='/panfs/roc/groups/4/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_wb_c='/panfs/jay/groups/2/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_template_nets='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';
       
    case '/tmp' %mesabi
        path{1}='/home/faird/shared/code/external/utilities/gifti-1.6';
        path{2}='/home/faird/shared/code/internal/utilities/Matlab_CIFTI';
        path{3}='/home/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/';
        path{4}='/home/faird/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';
        %templates:
        path{5}='/home/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii';
        path{6}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned.dscalar.nii';
        path{7}='/panfs/jay/groups/6/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii';
        path{8}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii'; 
        path{9}='/home/faird/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';        
        path{10}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dtseries.nii';
        path{11}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dtseries.nii';
        path{12}='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii';      
        path{13}='/home/faird/shared/code/internal/utilities/figure_maker/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene';      
        path{14}='/home/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_scalar_MSI.scene'; % avoid using the following, since it is intended for a label file instead of dscalar file. '/panfs/jay/groups/6/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_label_MSI.scene';      
        path{15}='/home/faird/shared/code/internal/utilities/figure_maker/make_dscalar_pics_v9.4.sh';
        path{16}='/home/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
        
        path_wb_c='/home/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; % workbench command path
        path_template_nets='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';
       
    case 'PCWIN64'
        path{1}='P:\code\external\utilities\gifti-1.6';
        path{2}='P:\code\external\utilities\Matlab_CIFTI';   
        path_wb_c='C:\Program Files\workbench-windows64-v1.1.1\workbench\bin_windows64\wb_command'; %path to wb_command
    otherwise
        error('Cannot parse directory structure to add dependencies automatically.')
end

settings.path=path;
settings.path_template_nets=path_template_nets;
settings.path_wb_c=path_wb_c;
