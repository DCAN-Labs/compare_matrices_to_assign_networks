function [all_nets_mat, all_nets_vec] = group_subcortical_network_proportions(dscalar_conc,outname,surface_only)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti/'))

%reference_dscalar = '/mnt/rose/shared/projects/ADHD_comm_det/ADHD_templmatch/rushmore_average/all_trio_and_prisma_TM_cleaned_5minutes_Control_avg.dscalar.nii';
%reference_dscalar = '/home/exacloud/lustre1/fnl_lab/projects/ADHD_comm_det/ADHD_templmatch/rushmore_average/all_trio_and_prisma_TM_cleaned_5minutes_Control_avg.dscalar.nii';
reference_dscalar = '/home/rando149/shared/projects/ABCD_net_template_matching/ABCD_GROUP_AVERAGES/template_matching/ABCD_group1_AVG_TM_Zscored_recolored.dscalar.nii'

%reference_dscalar = '/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.dscalar.nii';
%reference_dscalar = 'rushmore_average/all_trio_and_prisma_TM_cleaned_5minutes_Control_avg.dscalar.nii';


%Oregon_ADHD
%subject_list = importdata('/mnt/rose/shared/projects/ADHD_comm_det/ADHD_templmatch/rushmore_average/all_trio_and_prisma_TM_cleaned_5minutes_BOTH_rushpaths_dscalars.conc');
%NiGG Twins
%subject_list = importdata('/mnt/rose/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/template_matching_dscalars/template_matching_cleaned_dscalar.conc');
%ABCD Group1
%subject_list = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_group1_cleaned_TM_dscalars_trimmed.conc');
subject_list = importdata('/home/rando149/shared/projects/ABCD_net_template_matching/surfacearea/ABCD_templ_matched_scalars_group1_10_min_MSI_dscalars_trimmed.conc')

subject_list = importdata(dscalar_conc);
run_locally =0;

if run_locally ==1
%Some hardcodes:
wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
else

this_code = which('template_matching_RH');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

disp('Attempting to add neccesaary paths and functions.')
warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end
addpath(genpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti')) % remove non-working gifti path included with MSCcodebase
%rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
addpath(genpath('/mnt/max/shared/code/internal/utilities/plotting-tools'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command
end


if surface_only ==0
    all_nets_mat = zeros(size(subject_list,1), 21, 14);
    all_nets_vec = zeros(size(subject_list,1), 21*14);
else
    all_nets_mat = zeros(size(subject_list,1), 2, 14); %cortex only (i.e. only get proportions for left and right.)
    all_nets_vec = zeros(size(subject_list,1), 2*14);
end

for sub = 1: size(subject_list,1)
    disp(sub);
    %try
    %network_alluvial(DscalarC,DscalarD,plot_subcortex,output_name,wb_command)
    [~, structure_percentD] = network_alluvial(reference_dscalar,subject_list{sub},1,outname,wb_command);
    %catch
        structure_percentD =nan(21,14); % if there's an erro with the dscalar, fill with nans.
    %end
    all_nets_mat(sub,:,:) = structure_percentD;
    structure_percentD_trans = structure_percentD';
    structure_percentD_vec = structure_percentD_trans(:);
    structure_percentD_vec = structure_percentD_vec';
    all_nets_vec(sub,:) = structure_percentD_vec;
end
disp('Saving...')
%save('all_trio_and_prisma_TM_cleaned_5minutes_BOTH_dscalars_subcortical_proportion.mat','all_nets_mat','all_nets_vec')
save([outname '_subcortical_network_proportion.mat'],'all_nets_mat','all_nets_vec')

% [coeff,score,latent,tsquared,explained,mu] =pca(all_nets_vec);
% plot(1:size(explained,1),cumsum(100*explained(:,1)),'-bo');
% for i = 1:12
% figure()
% first2 = score(:,i:i+1);
% H = biplot(first2,'VarLabels',twin_str);
% end
% [XL,Yl,XS,YS,beta,PCTVAR] = plsregress(all_nets_vec,score,size(score,2));
% plot(1:size(score,2),cumsum(100*PCTVAR(2,:)),'-bo');
% imagesc(score)
% text(score(:,20),score(:,21),twin_str)
% twin_names_a = [1:2:25]';
% twin_names_b = [2:2:26]';
% 
% twin_str_a = num2str(twin_names_a');
% twin_str_b = num2str(twin_names_b');
% 
% figure()
% for i = 1:size(score,2)-1
% subplot(5,5,i)
% scatter(score(:,i),score(:,i+1))
% 
% text(score(1:2:(size(score,2)),i),score(1:2:size(score,2),i+1),twin_str_a,'Color','b')
% text(score(2:2:(size(score,2)+1),i),score(2:2:(size(score,2)+1),i+1),twin_str_b,'Color','r')
% end
disp('Done collecting network proportions.')
end

