%%%%%addpath(genpath('/projects/standard/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6'));
%%%%%addpath(genpath('/projects/standard/faird/shared/code/external/utilities/xmltree-2.0'));

%%%%%Example below
%%%%%mcc -v -m -d /mnt/max/shared/code/internal/utilities/corr_pt_dt -o cifti_conn_matrix_to_corr_pt_dt_exaversion /mnt/max/shared/code/internal/utilities/corr_pt_dt/cifti_conn_matrix_to_corr_pt_dt_exaversion.m -a /mnt/max/shared/code/external/utilities/Matlab_CIFTI/ -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m -a /mnt/max/shared/code/internal/utilities/corr_pt_dt/corr_pt_dt_exaversion.m

%%%%%system('/panfs/roc/msisoft/matlab/R2019a/bin/mcc -v -m -d /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -o template_matching_RH /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH.m -a /projects/standard/miran045/moral453/Desktop/container_template_matching/support_code/*.m -a /projects/standard/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m -a /projects/standard/faird/shared/code/external/utilities/xmltree-2.0/ -a /projects/standard/faird/shared/code/internal/utilities/plotting-tools/custom_hist -a /projects/standard/faird/shared/code/internal/utilities/Zscore_dconn -a /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -a /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files')

%%%%% mcc -v -m -R -singleCompThread -d /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -o template_matching_RH2 /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH_noaddpath4.m -a /projects/standard/miran045/moral453/Desktop/container_template_matching/support_code/*.m -a /projects/standard/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m -a /projects/standard/faird/shared/code/internal/utilities/plotting-tools/custom_hist -a /projects/standard/faird/shared/code/internal/utilities/Zscore_dconn -a /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -a /projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files -a /projects/standard/faird/shared/code/internal/utilities/plotting-tools -a /projects/standard/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6 -a /projects/standard/faird/shared/code/external/utilities/xmltree-2.0 % cristianMC added this 06/02/22

out_folder='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks'
out_name='template_matching_RH2'
input_function='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH_noaddpath5.m'

command='mcc -v -m -R singleCompThread';
command=strcat(command,{' '},{'-d '},out_folder);
command=strcat(command,{' '},{'-o '},out_name);
command=strcat(command,{' '},input_function);

in1='/projects/standard/miran045/moral453/Desktop/container_template_matching/support_code/*.m';
in2='/projects/standard/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m';
in3='/projects/standard/faird/shared/code/internal/utilities/plotting-tools/custom_hist';
in4='/projects/standard/faird/shared/code/internal/utilities/Zscore_dconn';
in5='/projects/standard/miran045/moral453/Desktop/container_template_matching/support_code/';
in6='/projects/standard/faird/shared/code/internal/utilities/plotting-tools';
in7='/projects/standard/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6';
in8='/projects/standard/faird/shared/code/external/utilities/xmltree-2.0';
in9='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files';
in10='/projects/standard/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/';


command=strcat(command,{' '},{'-a '},in1);
command=strcat(command,{' '},{'-a '},in2);
command=strcat(command,{' '},{'-a '},in3);
command=strcat(command,{' '},{'-a '},in4);
command=strcat(command,{' '},{'-a '},in5);
%%%%command=strcat(command,{' '},{'-a '},in6);
command=strcat(command,{' '},{'-a '},in7);
command=strcat(command,{' '},{'-a '},in8);
command=strcat(command,{' '},{'-a '},in9);
command=strcat(command,{' '},{'-a '},in10);

disp(command)

eval(char(command))


