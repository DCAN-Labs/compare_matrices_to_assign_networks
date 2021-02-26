function  make_overlap_distribution_figure
%-- 02/02/2020 01:46:26 PM --%
% This function loads in a matrix of network (16) x eta sqared value to
% network (91282).
load('Example_ABCD_eta_to_template_vox.mat')


this_code = which('template_matching_RH');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);
disp('Attempting to add neccesaary paths and functions.')
warning('off') %supress addpath warnings to nonfolders.
for i=2:np
addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};


histogram(eta_to_template_vox(:,1))
h = histogram(eta_to_template_vox(:,1));
load('PowerColorMap.mat')
h.FaceColor = [1 0 0];
h.FaceColor = 'r';
h.EdgeColor = 'r';
h.FaceAlpha = 1;
set(gcf,'color','w');
set(gca,'FontSize',20)



DMN_etas = eta_to_template_vox(:,1);

cii = ciftiopen(settings.path{8},wb_command);
cii.cdata = DMN_etas;
outputname  = 'Example_DMN_eta_to_template_per_voxel.dscalar.nii';
ciftisave(cii,outputname,wb_command)

cmd = (['/home/exacloud/lustre1/fnl_lab/code/internal/utilities/make_dscalar_pics/make_dscalar_pics_v8.sh ' pwd filesep 'Example_DMN_eta_to_template_per_voxel.dscalar.nii'  ' Example_DMN_eta_to_template_per_voxel ' pwd ' TRUE 0 0.75 videen_style FALSE none none none TRUE TRUE']);
disp(cmd)
system(cmd)

cii_dtseries = ciftiopen('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_templ_matched_scalars_group1_10_min/sub-NDARINVZR16R6Y3_ses-baselineYear1Arm1_task-rest_bold_timeseries_template_matched_Zscored_overlap_smooth_then_derivative_recolored.dtseries.nii',wb_command);
dtseries = cii_dtseries.cdata;

cii.cdata = dtseries(:,1);
outputname  = 'Example_DMN_labeled.dscalar.nii';
ciftisave(cii,outputname,wb_command)

cmd = (['/home/exacloud/lustre1/fnl_lab/code/internal/utilities/make_dscalar_pics/make_dscalar_pics_v8.sh ' pwd filesep 'Example_DMN_labeled.dscalar.nii'  ' Example_DMN_labeled ' pwd ' TRUE 1 18 power_surf TRUE 0.5 2 none TRUE TRUE']);
disp(cmd)
system(cmd)

end

