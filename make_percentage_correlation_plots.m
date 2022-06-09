function [all_rho, all_pval, all_nonmatched_num, all_nonmatched_percent] = make_percentage_correlation_plots(group1_conc_file,group2_conc_file,output_path,outputfilename_summary_name)

%group1_conc = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/group_pics/probability_pics/overlap_grp1.conc');
%group2_conc = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/group_pics/probability_pics/overlap_grp2.conc');


%add paths
this_code = which('make_percentage_correlation_plots');
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
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')
wb_command=settings.path_wb_c; %path to wb_command
group1_conc = importdata(group1_conc_file);
group2_conc = importdata(group2_conc_file);

%[~,outputfilename_summary_name] = fileparts(group1_conc_file);
for i=1:size(group1_conc)
 [~,B_temp,~] = fileparts(group1_conc{i});
 [~, outputfilename,~] = fileparts(B_temp);
%load data
%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_10min_GRP1_singlenet_percentage_n2988_DMN_network_percentage.dscalar.nii',wb_command);
%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_10min_GRP1_singlenet_probability_n2988_DMN_network.dscalar.nii',wb_command);

%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_GRP1_overlap_percentage_DMN_network_percentage.dscalar.nii',wb_command);
%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/group_pics/probability_pics/ABCD_GRP1_overlap_probability_DMN_network.dscalar.nii',wb_command);
cii = ciftiopen(group1_conc{i},wb_command);
dscalar_GRP1 = cii.cdata;

%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_10min_GRP2_singlenet_percentage_n3084_DMN_network_percentage.dscalar.nii',wb_command);
%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_10min_GRP2_singlenet_probability_n3084_DMN_network.dscalar.nii',wb_command);

%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/ABCD_percentage_maps/ABCD_GRP2_overlap_percentage_DMN_network_percentage.dscalar.nii',wb_command);
%cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/group_pics/probability_pics/ABCD_GRP2_overlap_probability_DMN_network.dscalar.nii',wb_command);

cii = ciftiopen(group2_conc{i},wb_command);
dscalar_GRP2 = cii.cdata;

% check here to make sure that the dscalars that are loaded in are the same
% length.  If one has a surface and one doesn't it will only compare
% surfaces.
if size(dscalar_GRP1,1) ~= size(dscalar_GRP2,1)
    if size(dscalar_GRP1,1) ==59412
        dscalar_GRP2 =dscalar_GRP2(1:59412,1);
    end
    if size(dscalar_GRP2,1) ==59412
        dscalar_GRP1 =dscalar_GRP1(1:59412,1);
    end
end

[rho_wzero,pval_wzero] = corr(dscalar_GRP1,dscalar_GRP2);

%get indices of nonzero elements.
log_nonZgrp1 = (dscalar_GRP1 ~=0);
log_nonZgrp2 = (dscalar_GRP2 ~=0);

both_nonzero = log_nonZgrp1 & log_nonZgrp2;

x = (dscalar_GRP1(both_nonzero));
y=(dscalar_GRP2(both_nonzero));

c = ksdensity([x,y], [x,y]);

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = [output_path filesep];
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 20;

% scaling
%fig.Units               = 'centimeters';
%fig.Position(3)         = opts.width;
%fig.Position(4)         = opts.height;


figure();scatter(x, y, 10,c,'filled');


set(gcf,'color','white')
set(gca,'FontSize',20)
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
%xlim([0 1]);ylim([0 1]);
fig.PaperPositionMode   = 'auto';
title('Correlation without Zeros')
[rho,pval] = corr(x,y);

log_Zerogrp1 = (dscalar_GRP1 ==0);
log_Zerogrp2 = (dscalar_GRP2 ==0);
disp(['The correlation between vectors is: ' num2str(rho) ', p=' num2str(pval) ])
disp(['The correlation between vectors including zero elements is: ' num2str(rho_wzero) ', p=' num2str(pval_wzero) ])

non_match_zero = xor(log_Zerogrp1,log_Zerogrp2);
disp(['Number of greyordinates with mismatched zeros is: ' num2str(sum(non_match_zero)) ' (or ' num2str((sum(non_match_zero))/size(both_nonzero,1)) '% of the greyordinates).'])

disp(['Saving image: ' opts.saveFolder outputfilename])
print([opts.saveFolder outputfilename], '-dpng', '-r600')


all_rho(i) = rho;
all_pval(i) = pval;
all_rho_wzero(i) = rho_wzero;
all_pval_wzero(i) = pval_wzero;

all_nonmatched_num(i) = sum(non_match_zero);
all_nonmatched_percent(i) = (sum(non_match_zero))/size(both_nonzero,1);
% nonZgrp1 = find(DMN_GRP1);
% nonZgrp1 = nnz(DMN_GRP1);
% nonZgrp1 = find(DMN_GRP1);
% nonZgrp1 = DMN_GRP1(find(DMN_GRP1));
% nonZgrp1 = DMN_GRP2(find(DMN_GRP2));
end
disp('Saving .mat file with summary stats saving...')
save([output_path filesep outputfilename_summary_name '_correlation_.mat'],'all_rho','all_pval','all_nonmatched_num','all_nonmatched_percent','all_rho_wzero','all_pval_wzero')
disp('Done running code for all files in conc.')

end