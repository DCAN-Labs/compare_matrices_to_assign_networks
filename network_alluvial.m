function [structure_percentC, structure_percentD] = network_alluvial(DscalarC,DscalarD,wb_command)

%OPTIONS
close all
check_surface =1;
check_subcortical =1;
sort_matrix =0;
plot_subcort =1;
%subcort_plot_type = 'alluvial';
subcort_plot_type = 'donut';
save_subcortical_percentages =1;
run_locally =0;
skip_plotting =1;
%left_labels = {'AUD', 'CO','DAN','DMN','FP','MTL','PMN','PON','SAL', 'SMD','SML' , 'Tpole', 'VAN','VIS'};
all_labels = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};

%right_labels = {'AUD', 'CO','DAN','DMN','FP','MTL','PMN','PON','SAL', 'SMD','SML' , 'Tpole', 'VAN','VIS'};
%right_labels = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};
possible_net_nums = [1 2 3  5  7 8 9 10 11 12 13 14 15 16];

% if run_locally ==1
% %Some hardcodes:
% wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
% addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
% else
% 
% this_code = which('template_matching_RH');
% [code_dir,~] = fileparts(this_code);
% support_folder=[code_dir '/support_files']; %find support files in the code directory.
% addpath(genpath(support_folder));
% settings=settings_comparematrices;%
% np=size(settings.path,2);
% 
% disp('Attempting to add neccesaary paths and functions.')
% warning('off') %supress addpath warnings to nonfolders.
% for i=2:np
%     addpath(genpath(settings.path{i}));
% end
% addpath(genpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti')) % remove non-working gifti path included with MSCcodebase
% %rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
% %rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
% addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
% addpath(genpath('/mnt/max/shared/code/internal/utilities/plotting-tools'));
% warning('on')
% wb_command=settings.path_wb_c; %path to wb_command
% end

%% Load data
%load('C:\Users\hermosir\Desktop\Test_data\MSC02a_to_ADHD315_template_MSC02b1_method_template_matching.mat')
%eta_subject_index_1min = eta_subject_index;
%A = eta_subject_index_1min;
%load('C:\Users\hermosir\Desktop\Test_data\MSC02a_to_ADHD315_template_MSC02b2_method_template_matching.mat')
%eta_subject_index_2min = eta_subject_index;
%B = eta_subject_index_2min;
%load('C:\Users\hermosir\Desktop\Test_data\MSC02a_to_ADHD315_template_MSC02b3_method_template_matching.mat')
%eta_subject_index_3min = eta_subject_index;
%C = eta_subject_index_3min;
%Ccii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\all_trio_and_prisma_TM_cleaned_5minutes_Control_avg.dscalar.nii',wb_command);
Ccii = ciftiopen(DscalarC,wb_command);

C = Ccii.cdata;
if check_subcortical ==1
    %C_info = ft_read_cifti_mod('C:\Users\hermosir\Documents\test_ciftis\all_trio_and_prisma_TM_cleaned_5minutes_Control_avg.dscalar.nii');
        C_info = ft_read_cifti_mod(DscalarC);

    greybrainstructs_C = C_info.brainstructure(find(C_info.brainstructure >0));
    surfacegreys_C = size(find(C_info.brainstructure <3 & C_info.brainstructure >0));
    C_surf = C(1:surfacegreys_C,1);
    C_sub = C(surfacegreys_C+1:end,1);
    for i = 1:size(C_info.brainstructurelabel,2)
        thisstructsgreys = find(greybrainstructs_C ==i);
        subcort_assings_by_structC{i,1} = C(thisstructsgreys);
    end
end
%load('C:\Users\hermosir\Desktop\Test_data\MSC02a_to_ADHD315_template_MSC02b4_method_template_matching.mat')
%eta_subject_index_4min = eta_subject_index;
%Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\all_trio_and_prisma_TM_cleaned_5minutes_ADHD_avg.dscalar.nii',wb_command);
Dcii = ciftiopen(DscalarD,wb_command);
%Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\sub-33015a_task-rest_DCANBOLDProc_v4.0.0_Atlas_template_matched_Zscored.dscalar.nii',wb_command);
%Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\sub-100501_ses-20100430_avg_number_of_networks.dscalar.nii',wb_command);

D = Dcii.cdata;
if check_subcortical ==1
    %D_info = ft_read_cifti_mod('C:\Users\hermosir\Documents\test_ciftis\all_trio_and_prisma_TM_cleaned_5minutes_ADHD_avg.dscalar.nii');
        D_info = ft_read_cifti_mod(DscalarD);
    greybrainstructs_D = D_info.brainstructure(find(D_info.brainstructure >0));
    surfacegreys_D = size(find(D_info.brainstructure <3 & D_info.brainstructure >0));
    D_surf = D(1:surfacegreys_D,1);
    D_sub = D(surfacegreys_D+1:end,1);
    for i = 1:size(D_info.brainstructurelabel,2)
        thisstructsgreys = find(greybrainstructs_D ==i);
        subcort_assings_by_structD{i,1} = D(thisstructsgreys);
    end
end
%test_missing_nets=1;
if check_surface ==1
    %D(find(D==10 | D==8 )) =5; % set all 10s to 1.
    D=D_surf;C=C_surf;
end
[alluvial_matrix, num_greys_in_all_nets_a, num_greys_in_all_nets_b, unique_nets_C, unique_nets_D]= build_alluvial_matrix(C, D);

% alluvial_matrix(4,:) = [];
% alluvial_matrix(:,4) = [];
% alluvial_matrix(5,:) = [];
% alluvial_matrix(:,5) = [];

isanet_in_current_setC = ismember(possible_net_nums,unique_nets_C); % used to pull a subset of colors from the powercolor set.
isanet_in_current_setD = ismember(possible_net_nums,unique_nets_D);

left_labels = all_labels(isanet_in_current_setC);
right_labels = all_labels(isanet_in_current_setD);

%h= alluvialflow_change_colors(alluvial_matrix, left_labels, right_labels, "Network Change",0,[],[], isanet_in_current_setC, isanet_in_current_setD); % unsorted

%    for i = this_nets_idices
%        if eta_subject_index_2min(this_nets_idices(i)) ~= i
%        else
%       end
num_greys_in_nets_a = sum(alluvial_matrix,2);
[E,Eidx] = sort(num_greys_in_nets_a,'descend');
num_greys_in_nets_b = sum(alluvial_matrix,1);
[F,Fidx] = sort(num_greys_in_nets_b,'descend');

if skip_plotting ==0
    if sort_matrix ==0 %
        h= alluvialflow_change_colors(alluvial_matrix, left_labels, right_labels, string('Network Change'),0,[],[], isanet_in_current_setC, isanet_in_current_setD); % unsorted
    else
        %disp(' sort matrix option is selected. NOTE: network colors need fixing.')
        %[E,Eidx] = sort(num_greys_in_all_nets_a,'descend');
        figure()
        set(gcf,'color','w');
        sorted_alluvial_a = alluvial_matrix(:,Fidx);
        sorted_alluvial_b = sorted_alluvial_a(Eidx,:);
        sorted_labels_a = left_labels(Eidx);
        sorted_labels_b = right_labels(Fidx);
        h = alluvialflow_change_colors(sorted_alluvial_b, sorted_labels_a, sorted_labels_b, string('Network Change'),1,Eidx, Fidx, isanet_in_current_setC, isanet_in_current_setD);
    end
end

if plot_subcort ==1
    if sort_matrix ==0 %      
        figure()       
        set(gcf,'color','w');
    end
    structure_percentC = zeros(size(subcort_assings_by_structC,1),size(possible_net_nums,2));
    structure_percentD = zeros(size(subcort_assings_by_structD,1),size(possible_net_nums,2));
    
    for k = 1: size(subcort_assings_by_structD,1)
        greysforstructC = subcort_assings_by_structC{k};
        greysforstructD = subcort_assings_by_structD{k};
        
        if save_subcortical_percentages ==1
            net_holder = 1;
            for x = possible_net_nums
                found_stuffC = find(greysforstructC == x);
                found_stuffD = find(greysforstructD == x);
                if isempty(found_stuffC) ==1
                    found_stuffC =0;
                else
                    structure_percentC(k,net_holder) = size(found_stuffC,1)/size(greysforstructC,1);
                end
                if isempty(found_stuffD) ==1
                    found_stuffD =0;
                else
                    structure_percentD(k,net_holder) = size(found_stuffD,1)/size(greysforstructD,1);
                end
                net_holder = net_holder +1;
            end
        else
        end
        
        if k>2
            [alluvial_matrix, num_greys_in_all_nets_a, num_greys_in_all_nets_b, unique_nets_C, unique_nets_D]= build_alluvial_matrix(greysforstructC, greysforstructD);
            isanet_in_current_setC = ismember(possible_net_nums,unique_nets_C); % used to pull a subset of colors from the powercolor set.
            isanet_in_current_setD = ismember(possible_net_nums,unique_nets_D);
            left_labels = all_labels(isanet_in_current_setC);
            right_labels = all_labels(isanet_in_current_setD);
            structure_name = C_info.brainstructurelabel{1,k};
            
            
            
            if sort_matrix ==0 %
                if skip_plotting ==0
                    subplot(4,5,k-2)
                end
                %set(gca,'XColor', 'none','YColor','none')
                switch subcort_plot_type
                    case 'alluvial'
                        h= alluvialflow_change_colors(alluvial_matrix, left_labels, right_labels, structure_name,0,[],[], isanet_in_current_setC, isanet_in_current_setD); % unsorted
                    case 'donut'
                        num_greys_in_nets_a = sum(alluvial_matrix,2);
                        %                        [E,Eidx] = sort(num_greys_in_nets_a,'descend');
                        num_greys_in_nets_b = sum(alluvial_matrix,1);
                        num_greys_in_nets_b = num_greys_in_nets_b';
                        
                        %                        [F,Fidx] = sort(num_greys_in_nets_b,'descend');
                        dough{1,1} = num_greys_in_nets_a';
                        dough{2,1} = num_greys_in_nets_b';
                        network_donut(dough,left_labels, right_labels, structure_name,0,Eidx, Fidx, isanet_in_current_setC, isanet_in_current_setD,skip_plotting);
                        
                end
            else
                %disp(' sort matrix option is selected. NOTE: network colors need fixing.')
                %[E,Eidx] = sort(num_greys_in_all_nets_a,'descend');
                num_greys_in_nets_a = sum(alluvial_matrix,2);
                [E,Eidx] = sort(num_greys_in_nets_a,'descend');
                num_greys_in_nets_b = sum(alluvial_matrix,1);
                num_greys_in_nets_b = num_greys_in_nets_b';
                [F,Fidx] = sort(num_greys_in_nets_b,'descend');
                sorted_alluvial_a = alluvial_matrix(:,Fidx);
                sorted_alluvial_b = sorted_alluvial_a(Eidx,:);
                sorted_labels_a = left_labels(Eidx);
                sorted_labels_b = right_labels(Fidx);
                if skip_plotting ==0
                    
                    subplot(4,5,k-2)
                end
                switch subcort_plot_type
                    case 'alluvial'
                        h = alluvialflow_change_colors(sorted_alluvial_b, sorted_labels_a, sorted_labels_b, structure_name,1,Eidx, Fidx, isanet_in_current_setC, isanet_in_current_setD);
                    case 'donut'
                        dough{1,1} = num_greys_in_nets_a';
                        dough{2,1} = num_greys_in_nets_b';
                        network_donut(dough,left_labels, right_labels, structure_name,1,Eidx, Fidx, isanet_in_current_setC, isanet_in_current_setD,skip_plotting);  
                end
            end
        end %ifk>2
    end
end
disp('Done plotting alluvial diagrams')
% if overlap ==0
%     left_labels = {'DAN'};
%     %DAN_split = alluvial_matrix(3,:);
%     DAN_split = net_empty;
%
% end
% figure()
% g= alluvialflow(DAN_split, 'DAN', right_labels, "Network Change");

end

function [alluvial_matrix, num_greys_in_all_nets_a, num_greys_in_all_nets_b, unique_nets_C, unique_nets_D]= build_alluvial_matrix(vecC, vecD)
unique_nets_C =unique(vecC);
num_rows = size(unique_nets_C,1);
unique_nets_D = unique(vecD);
num_cols = size(unique_nets_D,1);
alluvial_matrix = zeros(num_rows,num_cols); % preallocate alluvial arrary (for template-matching usually a 14 x 14 matrix).

N=1; %make counters for the array
for n = unique_nets_C' % add 2 here to make the alluvial matrix 16 x 16 instead of 14 x 14.
    C_net_indices = find(vecC == n) ;
    num_greys_in_all_nets_a(N) = size(C_net_indices,1); % save the size of each network for later.
    M=1;
    for m = unique_nets_D'
        D_nets_indices = find(vecD == m) ;
        num_greys_in_all_nets_b(M) = size(D_nets_indices,1);
        alluvial_matrix(N,M) = sum((vecD(C_net_indices)) == m); % save the size of each network for later.
        M=M+1;
    end
    N=N+1;
end

% for n = 1: num_rows+2 % add 2 here to make the alluvial matrix 16 x 16 instead of 14 x 14.
%     this_nets_indices = find(vecC == n) ;
%     num_greys_in_all_nets_a(n) = size(this_nets_indices,1); % save the size of each network for later.
%     for m = 1:num_cols+2
%         other_nets_indices = find(vecD == m) ;
%         num_greys_in_all_nets_b(m) = size(other_nets_indices,1);
%         alluvial_matrix(n,m) = sum((vecD(this_nets_indices)) == m); % save the size of each network for later.
%     end
% end

%alluvial_matrix =alluvial_matrix';
end