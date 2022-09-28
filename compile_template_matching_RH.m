%%%%%addpath(genpath('/home/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6'));
%%%%%addpath(genpath('/home/faird/shared/code/external/utilities/xmltree-2.0'));

%%%%%Example below
%%%%%mcc -v -m -d /mnt/max/shared/code/internal/utilities/corr_pt_dt -o cifti_conn_matrix_to_corr_pt_dt_exaversion /mnt/max/shared/code/internal/utilities/corr_pt_dt/cifti_conn_matrix_to_corr_pt_dt_exaversion.m -a /mnt/max/shared/code/external/utilities/Matlab_CIFTI/ -a /mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/isthisanoutlier.m -a /mnt/max/shared/code/internal/utilities/corr_pt_dt/corr_pt_dt_exaversion.m

%%%%%system('/panfs/roc/msisoft/matlab/R2019a/bin/mcc -v -m -d /home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -o template_matching_RH /home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH.m -a /home/miran045/moral453/Desktop/container_template_matching/support_code/*.m -a /panfs/roc/groups/8/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m -a /home/faird/shared/code/external/utilities/xmltree-2.0/ -a /home/faird/shared/code/internal/utilities/plotting-tools/custom_hist -a /home/faird/shared/code/internal/utilities/Zscore_dconn -a /home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -a /home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files')

%%%%% mcc -v -m -R -singleCompThread -d /panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -o template_matching_RH2 /panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH_noaddpath4.m -a /panfs/roc/groups/4/miran045/moral453/Desktop/container_template_matching/support_code/*.m -a /panfs/roc/groups/8/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m -a /panfs/roc/groups/8/faird/shared/code/internal/utilities/plotting-tools/custom_hist -a /panfs/roc/groups/8/faird/shared/code/internal/utilities/Zscore_dconn -a /panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks -a /panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files -a /panfs/roc/groups/8/faird/shared/code/internal/utilities/plotting-tools -a /panfs/roc/groups/8/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6 -a /panfs/roc/groups/8/faird/shared/code/external/utilities/xmltree-2.0 % cristianMC added this 06/02/22

out_folder='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks'
out_name='template_matching_RH2'
input_function='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/template_matching_RH_noaddpath5.m'

command='mcc -v -m -R singleCompThread';
command=strcat(command,{' '},{'-d '},out_folder);
command=strcat(command,{' '},{'-o '},out_name);
command=strcat(command,{' '},input_function);

in1='/home/miran045/moral453/Desktop/container_template_matching/support_code/*.m';
in2='/home/faird/shared/code/internal/utilities/cifti_connectivity/src/*.m';
in3='/home/faird/shared/code/internal/utilities/plotting-tools/custom_hist';
in4='/home/faird/shared/code/internal/utilities/Zscore_dconn';
in5='/home/miran045/moral453/Desktop/container_template_matching/support_code/';
in6='/home/faird/shared/code/internal/utilities/plotting-tools';
in7='/home/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6';
in8='/home/faird/shared/code/external/utilities/xmltree-2.0';
in9='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files';
in10='/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/';


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


