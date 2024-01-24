function [final_patch_path] = patch_match(subject_input_cifti_file,template_input_cifti_file,output_template_path,path_to_Lmidthicknessfile,path_to_Rmidthicknessfile, output_subject_path,output_file_name,distance_matrix_to_use)

%This function work trying to match up a subjects individualized clusters
%with the clusters observed in the group.

%R. Hermosillo 09/29/2021

% inputs are:
% "subject_input_cifti_file" =  the full path to the subject dscalar with their template matching networks that you
%                               are going to try to match against some template.
% "template_input_cifti_file" = this is template that you wish to match
%                               the cluster to.  ABCD has around 200 patches depending on how the
%                               clusters are defined using the parameteres specified below.  Use 'default
%                               to use the ABCD_consesus map.
% "output_template_path" =      path to where you want to write the template patches file. 
% "path_to_Lmidthicknessfile" = path to the LEFT midthickness file that you will use to find the clusters
% "path_to_Rmidthicknessfile" = path to the RIGHT midthickness file that you will use to find the clusters
% "output_subject_path" =       path to output the subject's files. Several 
%                               intermediate files are written (e.g.
%                               patches are identified per network, and
%                               the final files).
% "output_file_name" =          Some output file name to use (E.g. "SUBJECT1234")
% "distance_matrix_to_use" =    Path to a  geodesic+euclidean distance 
%                               matrix.  If you haven't created one before you can use the code found
%                               here: /home/faird/shared/code/internal/utilities/distance-matrix or 
%                                at https://gitlab.com/Fair_lab/distance-matrix.git
% 
% 
% An Example call to this functin would be as follows:
% patch_match('/path/to/my/TM_networks.dscalar.nii','default','/path/to/my/output_folder/,'/path/to/my/L_midthickness_file.surf.gii','/path/to/my/R_midthickness_file.surf.gii','the_best_subject_ever','/path/to/my/huge/distance_matrix_file.mat');

% HOW DOES THE CODE WORK?

%First both small islands are removed from both the template and the
%subject using a size exclusion or 80mm^3.  because combinations of patches
%can be used, don't consider clusters that are less than 10 grayordinates
%in this calculation.

%In the first step cluster files are created for each network.  Then the
%jaccard similiarty is calculated for each patch again sthe the patches in
%the template.  AFter that, combinations of patches are also tested against
%the template patches.  Most grayordinates are assigned at this step.

%In the second step For small islands Distance to unassigned group-patches of the same
%network.  This allows for patches that are very close to the template (but
%maybe missed in step 1) to be assigned.  The code works by calculating the
%distance between each grayordinate to the unassigned template patches.

%Lastly, The code the assigns the remaining small islands to the patches that have already been previously assigned based on distance.
%If the distance is too large, those islands will remain unassigned.



%% Start with some hardcoded parameters.
tic
%Hardcodes
%subject_input_cifti_file='C:\Users\hermosir\Documents\test_ciftis\sub-33005b_task-rest_DCANBOLDProc_v4.0.0_Atlas_template_matched_Zscored_recolored.dscalar.nii';
%template_input_cifti_file=('C:\Users\hermosir\Documents\repos\support_folder\91282_Greyordinates.dscalar.nii');
%template_input_cifti_file=('C:\Users\hermosir\Documents\test_ciftis\ABCD_GROUP_AVERAGES\ABCD_group1_AVG_TM_Zscored_recolored.dscalar.nii');
%template_input_cifti_file='C:\Users\hermosir\Documents\repos\support_folder\ABCD_GRP1_91282_Greyordinates_consensus_recolored.dscalar.nii';
switch template_input_cifti_file
    case 'default'
        template_input_cifti_file='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/ABCD_GRP1_91282_Greyordinates_consensus_recolored.dscalar.nii';
    otherwise
        disp(['Using template file: ' template_input_cifti_file]);
end
%output_template_path = 'C:\Users\hermosir\Documents\repos\support_folder\ABCD_GRP1_avg_network_patches';
%output_subject_path = 'C:\Users\hermosir\Documents\repos\support_folder\ABCD_GRP1_avg_network_patches';
%output_file_name = 'sub-33005b';
min_patch_size = 80; % previously 30, but since bold voxels are 2x2x2, only 4 voxels=32.  80 means that voxels must be at least 10 grayodrinates
min_num_of_grays = 10; % with above, if a cluster has less than 10 grayordiantes, don't count it as a unique patch to match (for either the template or the subject).
maximum_combination_of_nets = 4;
save_matched_dscalars =1;
min_dist=30;
keep_cortical_subcortical_seperation =1; %Set to 1 to set cortico-subcortical distance at 255mm (max).  If 0, distance matrix will use the eucliden distance from cortical to subcortical grayordiantes.

%distance_matrix_to_use = [support_folder filesep 'EUGEODistancematrix_XYZ_255interhem_unit8.mat'];
%distance_matrix_to_use = [support_folder filesep 'EUGEODistancematrix_XYZ_unit8.mat']

%     if exist('distances','var') ==1
%         disp('Distance matrix already loaded.')
%     else
%         tic
%         disp('loading distance matrix...')
%         if keep_cortical_subcortical_seperation ==1
%             disp('Note: Distance matrix has eucliean distances between the cortex and subcortex set to 255mm (max uint8).')
%             load([support_folder filesep 'EUGEODistancematrix_XYZ_255interhem_unit8.mat'],'distances');
%         else
%             disp('Note: Distance matrix uses eucliean distances between the cortex and subcortex.')
%             load([support_folder filesep 'EUGEODistancematrix_XYZ_unit8.mat'],'distances');
%         end
%         toc
%     end


%% Step 0: Add dependency paths
run_locally =0;
%add cifti paths
if run_locally ==1
    %Some hardcodes:
    wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
    support_folder='C:\Users\hermosir\Documents\repos\support_folder';
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
    warning('on')
    if exist('wb_command','var') ==1
    else
        wb_command=settings.path_wb_c; %path to wb_command
    end
end

%load power colors
%load('C:\Users\hermosir\Documents\repos\support_folder\Jet_wzerowhite_colormap.mat','mymap')

net_list = [1 2 3 5 7 8 9 10 11 12 13 14 15 16]; % hardcoded network assingments.
%all_labels = {'DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};

%% Step 1: Load template data and create cluster_files
% Let's get started
if exist([output_template_path filesep 'template_net_all_unique_patches.dscalar.nii'],'file')
    load([output_template_path filesep 'template_net_all_unique_patches.mat']);
else
    disp('Template data cannot be found.')
    %error('You probably do not want to remake the template')
    template_cifti_obj = ciftiopen(template_input_cifti_file,wb_command);
    template_nets = template_cifti_obj.cdata;
    
    disp('Creating_patch_file_from_template_dscalar');
    disp(['Minimum patch size = ' num2str(min_patch_size)])
    k=1;
    
    for net_num =net_list
        disp('Making patch dscalars for each network for template...')
        template_raw_seperated_cifti_file = [output_template_path filesep 'template_net_' num2str(net_num) '_raw.dscalar.nii'];
        patch_outputname_cifti_file = [output_template_path filesep 'template_net_' num2str(net_num) 'patches_size' num2str(min_patch_size) '.dscalar.nii'];
        
        patch_list{k,1} = patch_outputname_cifti_file;k=k+1;
        
        if exist(patch_outputname_cifti_file,'file') ~= 0
            disp('Template patch dscalars has already been made for this network.')
        else
            
            %net_num = 7;
            this_net = template_nets==net_num;
            this_net_double = double(this_net);
            template_cifti_obj.cdata=this_net_double;
            ciftisave(template_cifti_obj,template_raw_seperated_cifti_file,wb_command);
            
            cmd=[wb_command ' -cifti-find-clusters ' template_raw_seperated_cifti_file ' 0 ' num2str(min_patch_size) ' 0  ' num2str(min_patch_size) ' COLUMN ' patch_outputname_cifti_file ' -left-surface ' path_to_Lmidthicknessfile  ' -right-surface ' path_to_Rmidthicknessfile  ];
            %disp(cmd);
            system(cmd);
        end
    end
    
    %Step 2 Load patches
    disp('loading patch dscalars...')
    for k=1:size(patch_list,1)
        template_patch_obj = ciftiopen(patch_list{k},wb_command);
        patch_matrix_column = net_list(k);
        patch_matrix(:,patch_matrix_column) = template_patch_obj.cdata;
    end
    
    adjusted_patch_matrix = zeros(size(patch_matrix,1),size(patch_matrix,2));
    for j = 1: size(patch_matrix,2)
        if j ==1
            patches_this_net = nonzeros(unique(patch_matrix(:,j)));
            full_patches_list = patches_this_net;
            unique_patches_so_far = nonzeros(unique(full_patches_list));
            template_list_of_adjusted_patches_by_networks{j,1} = patches_this_net;
            template_list_of_patches_by_networks{j,1} = patches_this_net;
            template_num_patches(j,1) = length(nonzeros(unique(full_patches_list)));
            for k = 1:size(patches_this_net,1)
                this_patches_indices = find(patch_matrix(:,j)==patches_this_net(k));
                adjusted_patch_matrix(this_patches_indices,j) = patches_this_net(k);
            end
            
        else
            patches_this_net = nonzeros(unique(patch_matrix(:,j)));
            if isempty(unique_patches_so_far) ==1
            else
                for k = 1:size(patches_this_net,1)
                    this_patches_indices = find(patch_matrix(:,j)==patches_this_net(k));
                    adjusted_patch_matrix(this_patches_indices,j) = patches_this_net(k)+max(full_patches_list);
                end
                adjusted_patches_this_net = patches_this_net+max(full_patches_list);
                full_patches_list = [full_patches_list; adjusted_patches_this_net];
                template_list_of_adjusted_patches_by_networks{j,1} = adjusted_patches_this_net;
                template_list_of_patches_by_networks{j,1} = patches_this_net;
                template_num_patches(j,1) = length(nonzeros(unique(full_patches_list)));
                %template_list_of_patches_by_networks{j,1} = patches_this_net;
            end
        end
    end

template_patch_label_vector =sum(adjusted_patch_matrix,2);
template_cifti_obj.cdata=template_patch_label_vector;
ciftisave(template_cifti_obj,[output_template_path filesep 'template_net_all_unique_patches.dscalar.nii'],wb_command);

num_greys_left = nnz(template_patch_label_vector);
num_greys_removed = size(adjusted_patch_matrix,1)-nnz(template_patch_label_vector);
disp(['Number of grayordinates removed: ' num2str(num_greys_removed)]);
disp(['Number of grayordinates left: ' num2str(num_greys_left)]);
disp(['Number of patches at this threshold: ' num2str(max(full_patches_list))]);


end


%% Step 3 :Load subject data.
disp(['Loading subject...' subject_input_cifti_file]);
subject_cifti_obj = ciftiopen(subject_input_cifti_file,wb_command);
subject_nets = subject_cifti_obj.cdata;

disp('Creating patch file from subject dscalar');
disp(['Minimum patch size = ' num2str(min_patch_size)])
k=1;

for net_num =net_list
    disp('Making patch dscalars for each network for subject...')
    subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_raw.dscalar.nii'];
    patch_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) 'patches_size' num2str(min_patch_size) 'patches.dscalar.nii'];
    subject_patch_list{k,1} = patch_subject_outputname_cifti_file;k=k+1;
    
    if exist(patch_subject_outputname_cifti_file,'file') ~=0
        disp('Subject patch dscalars has already been made for this network.')
    else
        
        %net_num = 7;
        this_net = subject_nets==net_num;
        this_net_double = double(this_net);
        subject_cifti_obj.cdata=this_net_double;
        ciftisave(subject_cifti_obj,subject_outputname_cifti_file,wb_command);
        
        cmd=[wb_command ' -cifti-find-clusters ' subject_outputname_cifti_file ' 0 ' num2str(min_patch_size) ' 0  ' num2str(min_patch_size) ' COLUMN ' patch_subject_outputname_cifti_file ' -left-surface ' path_to_Lmidthicknessfile ' -right-surface ' path_to_Rmidthicknessfile ];
        %disp(cmd);
        system(cmd);
    end
end

%Step 2 Load patches
disp('loading patch dscalars...')
for k=1:size(subject_patch_list,1)
    subject_patch_obj = ciftiopen(subject_patch_list{k},wb_command);
    patch_matrix_column = net_list(k);
    subject_patch_matrix(:,patch_matrix_column) = subject_patch_obj.cdata;
end

subject_adjusted_patch_matrix = zeros(size(subject_patch_matrix,1),size(subject_patch_matrix,2));

for j = 1: size(subject_patch_matrix,2)
    if j ==1
        patches_this_net = nonzeros(unique(subject_patch_matrix(:,j)));
        subject_full_patches_list = patches_this_net;
        unique_patches_so_far = nonzeros(unique(subject_full_patches_list));
        subject_list_of_adjusted_patches_by_networks{j,1} = patches_this_net;
        for k = 1:size(patches_this_net,1)
            this_patches_indices = find(subject_patch_matrix(:,j)==patches_this_net(k));
            subject_adjusted_patch_matrix(this_patches_indices,j) = patches_this_net(k);
        end
        
    else
        patches_this_net = nonzeros(unique(subject_patch_matrix(:,j)));
        if isempty(unique_patches_so_far) ==1
        else
            for k = 1:size(patches_this_net,1)
                this_patches_indices = find(subject_patch_matrix(:,j)==patches_this_net(k));
                subject_adjusted_patch_matrix(this_patches_indices,j) = patches_this_net(k)+max(subject_full_patches_list);
            end
            adjusted_patches_this_net = patches_this_net+max(subject_full_patches_list);
            subject_list_of_adjusted_patches_by_networks{j,1} = adjusted_patches_this_net;
            subject_full_patches_list = [subject_full_patches_list; adjusted_patches_this_net];
        end
    end
end

subject_patch_label_vector =sum(subject_adjusted_patch_matrix,2);
subject_cifti_obj.cdata=subject_patch_label_vector;
ciftisave(subject_cifti_obj,[output_subject_path filesep output_file_name 'subject_net_all_unique_patches.dscalar.nii'],wb_command);

sub_num_greys_left = nnz(subject_patch_label_vector);
sub_num_greys_removed = size(subject_adjusted_patch_matrix,1)-nnz(subject_patch_label_vector);
disp(['Number of grayordinates removed: ' num2str(sub_num_greys_removed)]);
disp(['Number of grayordinates left: ' num2str(sub_num_greys_left)]);
disp(['Number of patches at this threshold: ' num2str(max(subject_full_patches_list))]);


%%step 4 Create logicals for each cluster separated by network.

for n =1: size(subject_list_of_adjusted_patches_by_networks,1)
    if n == 4 || n ==6
    else
        subject_adj_cluster_nums_for_this_net = subject_list_of_adjusted_patches_by_networks{n,1};
        template_adj_cluster_nums_for_this_net = template_list_of_adjusted_patches_by_networks{n,1};
        
        for p = 1:size(subject_adj_cluster_nums_for_this_net,1)
            subject_all_patches_this_network_separated(:,p)=(subject_patch_label_vector ==subject_adj_cluster_nums_for_this_net(p));
        end
        
        for p = 1:size(template_adj_cluster_nums_for_this_net,1)
            template_all_patches_this_network_separated(:,p)=(template_patch_label_vector ==template_adj_cluster_nums_for_this_net(p));
        end
        
        subject_all_clust_log_vectors{n}=subject_all_patches_this_network_separated;
        template_all_clust_log_vectors{n}=template_all_patches_this_network_separated;
        %clear these variables so that the cell sizes accurate array sizes within cells are preserved.
        clear subject_all_patches_this_network_separated template_all_patches_this_network_separated;
    end
end

%% step 5 Jaccard time.
jaccard_mat_all_nets=cell(size(subject_all_clust_log_vectors,2),1);
for n= 1:size(subject_all_clust_log_vectors,2)
    disp(['Matching patches for network ' num2str(n)])
    
    sub_num_clusters_in_this_net = size(subject_all_clust_log_vectors{n},2);
    template_num_clusters_in_this_net = size(template_all_clust_log_vectors{n},2);
    subject_this_nets_cluster_logicals = subject_all_clust_log_vectors{n}; %select this network
    template_this_nets_cluster_logicals = template_all_clust_log_vectors{n}; %select this network
    
    disp(['Matching patches for network ' num2str(n)])
    disp(['Number of network ' num2str(n) ' patches for template= '  num2str(template_num_clusters_in_this_net) ]);
    disp(['Number of network ' num2str(n) ' patches for subject = '  num2str(sub_num_clusters_in_this_net) ]);
    
    
    for p = 1:template_num_clusters_in_this_net
        patch2matchto = template_this_nets_cluster_logicals(:,p);
        disp(['Calculating match (Jaccard) for patch '  num2str(p) ' for network: ' num2str(n) ' ...' ]);
        D{n,1}=[]; k=1;
        
        this_patch_vec = 1:sub_num_clusters_in_this_net;
        for i = 1:maximum_combination_of_nets
            C = nchoosek(this_patch_vec,k); % return all possible combinations of select k nets from all networks.
            % put the vector of C into a matrix of zeros so that it can be concatenated onto the full list.
            Czeromat = zeros(size(C,1),size(this_patch_vec,2),1);
            Czeromat(:,1:size(C,2)) = C;
            D{n,1} = [D{n,1}; Czeromat];
            k=k+1;
        end
        if p==1
            jaccard_mat=zeros(size(D{n,1},1),size(template_this_nets_cluster_logicals,2));
        end
        for q = 1:size(D{n,1},1)
            net_indices_combo = nonzeros(D{n,1}(q,:));
            E=logical(sum(subject_this_nets_cluster_logicals(:,net_indices_combo),2));
            jaccard_mat(q,p)=jaccard(E,patch2matchto);
        end
    end
    jaccard_mat_all_nets{n}=jaccard_mat;
    
end

%Step 6 Set patch numbers with largest jaccard.

for n = 1: size(jaccard_mat_all_nets,1)
    if n == 4 || n ==6
    else
        jaccard_this_net = jaccard_mat_all_nets{n,1};
        
        %Before finding out which combos have the largest jaccard, check
        %that the the maximum values for jaccard are not all zero (which
        %indicates no overlap of any subject patches with the given template patch).
        %If any of the columnes in the jaccard matrix are all zero, then
        %matlab's max function will label the first element of the matrix
        %the maximum.  Here, we set them to Nan.
        
        for j=1:size(jaccard_this_net,2)
            jaccard_checksum = sum(jaccard_this_net,1);
            for k =1: size(jaccard_checksum,2)
                if jaccard_checksum(k) ==0
                    jaccard_mat_all_nets{n,1}(:,k)=nan;
                end
            end
        end
        
        %make copies of these for later. Because they will be modified to
        %remove patches from the pool. (Akin to sample without replacement)
        jaccard_mod = jaccard_mat_all_nets;
        possible_patch_values_template = 1:size(jaccard_mod{n,1},2);
        num_sub_clust= max(D{n,1}(:,1));
        assignment_mat{n,1} = zeros(num_sub_clust,1);
        
        for i = 1:size(jaccard_mod{n,1},2)
            maximum = max(max(jaccard_mod{n,1}));
            if maximum ==0
                disp('max is zero')
                break
            else
                [x,y]=find(jaccard_mod{n,1}==maximum);
                if size(x,1)>1
                    x=x(1);
                    y=y(1);
                end
                jaccard_mod{n,1}(x,:)=nan;
                jaccard_mod{n,1}(:,y)=nan;
                patches_to_exclude = nonzeros(D{n,1}(x,:));
                
                %assingn nets
                assignment_mat{n,1}(patches_to_exclude) = y;
                
                D_mod{n,1} = ismember(D{n,1}(:,:),patches_to_exclude); %find combos that have already been assigned.
                
                %exclude additional jaccard values if the network has already been
                %assigned.
                additional_combos_to_exclude = any(D_mod{n,1},2);
                jaccard_mod{n,1}(additional_combos_to_exclude,:)=nan;
                
                possible_patch_values_template(possible_patch_values_template==y) =[];
            end
        end
        %[max_values{n},patch_indices_this_net] =  max(jaccard_mat_all_nets{n,1},[],1,'includenan');
        %patch_indices_this_net(all(isnan(jaccard_mat_all_nets{n,1}),1))= NaN;
        
        %         for i=1:size(patch_indices_this_net,2)
        %             if isnan(patch_indices_this_net(i)) ==1
        %                 best_combos{n}(i,:) = NaN;
        %             else
        %                 best_combos{n}(i,:) = D{n,1}(patch_indices_this_net(i),:);
        %                 %best_combos{n}(i,:) = (D{n,1}(patch_indices_this_net,i));
        %
        %             end
        %         end
    end
end

for n=1:size(assignment_mat,1)
    unlabeled_patches(n,1) = size(assignment_mat{n,1},1)-nnz(assignment_mat{n,1});
    labeled_patches(n,1)= nnz(assignment_mat{n,1});
end


subject_patch_matched_matrix=zeros(size(subject_patch_matrix,1),size(subject_patch_matrix,2));


for i= 1:size(assignment_mat,1)
    for j=1:size(assignment_mat{i,1},1)
        orig_patch_dix=subject_patch_matrix(:,i)==j; % get the logical indices of each patch.
        subject_patch_matched_matrix(orig_patch_dix,i)=assignment_mat{i,1}(j);
    end
end

if save_matched_dscalars==1
    %k=1;
    for net_num =net_list
        disp('Saving patch-matched dscalars for each network for subject...')
        %subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '.dscalar.nii'];
        patch_matched_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched.dscalar.nii'];
        %subject_patch_matched_list{k,1} = patch_matched_subject_outputname_cifti_file;
        
        %if exist(patch_matched_subject_outputname_cifti_file,'file') ~=0
        %    disp('Subject patch dscalars has already been made for this network.')
        %else
        this_net_double = double(subject_patch_matched_matrix(:,net_num));
        subject_cifti_obj.cdata=this_net_double;
        ciftisave(subject_cifti_obj,patch_matched_subject_outputname_cifti_file,wb_command);
        %system(cmd);
        %end
        %k=k+1;
    end
else
end

total_unlabelled_patches = sum(unlabeled_patches);
total_labelled_patches = sum(labeled_patches);
patch_assignedgreys = any(subject_patch_matched_matrix,2);
number_of_grey_assigned_with_jaccard_only = sum(patch_assignedgreys);
disp(['The number of unlabelled patched with jaccard only is: ' num2str(total_unlabelled_patches) '/' num2str(max(subject_full_patches_list)) ]);
disp(['The number of   labelled patched with jaccard only is: ' num2str(total_labelled_patches) '/' num2str(max(full_patches_list)) ' template patches ']);
disp(['The number of greys assigned with jaccard only is: ' num2str(number_of_grey_assigned_with_jaccard_only) '/' num2str(size(patch_assignedgreys,1))]);
pause(4) % so that people can read
toc

disp('Continuing...')
disp('Next step: Assigning unassigned patches based on...')
disp(['Step 2/3 - Distance to unassigned group-patches of the same network using proximity: ' num2str(min_dist)])
pause(4)

for i= 1:size(assignment_mat,1)
    template_unmatched_clusters_poststep1{i,1} =  setdiff(template_list_of_patches_by_networks{i,1}(:),assignment_mat{i,1}(:));
end

for i=1:size(template_unmatched_clusters_poststep1,1)
    if isempty(template_unmatched_clusters_poststep1{i,1}) ~=1
        for j=1:size(template_unmatched_clusters_poststep1{i,1},1)
            template_unmatched_patch_indices{i,j} = find(patch_matrix(:,i)==template_unmatched_clusters_poststep1{i,1}(j));
        end
    end
end

i=1;k=1;
for i=1:size(assignment_mat,1)
    %if isempty(assignment_mat{i,1}) ~=1
    subject_unassigned_patches{i,1} = find(assignment_mat{i,1}==0);
    %adjusted_missing_nets{i,1}=subject_list_of_adjusted_patches_by_networks{i,1}(unassigned_patches);
    for k=1:size(subject_unassigned_patches{i,1},1)
        subject_adjusted_missing_nets_indices_by_net{i,k}=find(subject_patch_matrix(:,i)==subject_unassigned_patches{i,1}(k));
    end
    %end
end

if exist([output_subject_path filesep output_file_name '_templatepatchto' output_file_name 'patch_distance_mat_by_net.mat'],'file')==2
    disp('loading patch2patch distances...')
    load([output_subject_path filesep output_file_name '_templatepatchto' output_file_name 'patch_distance_mat_by_net.mat'],'templatepatchtosubjectpatch_distance_mat_by_net')
else
    
    if exist('distances','var') ==1
        disp('Distance matrix already loaded.')
    else
        tic
        disp('loading distance matrix...')
        if keep_cortical_subcortical_seperation ==1
            disp('Note: Distance matrix has eucliean distances between the cortex and subcortex set to 255mm (max uint8).')
            %load([support_folder filesep 'EUGEODistancematrix_XYZ_255interhem_unit8.mat'],'distances');
            load(distance_matrix_to_use,'distances');
        else
            disp('Note: Distance matrix uses eucliean distances between the cortex and subcortex.')
            %load([support_folder filesep 'EUGEODistancematrix_XYZ_unit8.mat'],'distances');
            load(distance_matrix_to_use,'distances');
        end
        toc
    end
    % calculate the distance from each unassigned template patch to unassigned subject
    % patches.
    for i=1:size(template_unmatched_patch_indices,1)
        if i ==1 && isempty(template_unmatched_patch_indices{i,j}) ==1
            templatepatchtosubjectpatch_distance_mat_by_net{i,1}=[];
        else
            num_unpatched_template_pathces=sum(~cellfun(@isempty,template_unmatched_patch_indices(i,:)));
            for j=1:num_unpatched_template_pathces
                template_lonely_patch_indices = template_unmatched_patch_indices{i,j};
                
                if isempty(template_lonely_patch_indices) ==1
                    templatepatchtosubjectpatch_distance_mat{j,1} = [];
                else
                    for k=1:size(subject_adjusted_missing_nets_indices_by_net,2)
                        subject_lonely_patch_indices = subject_adjusted_missing_nets_indices_by_net{i,k};
                        if isempty(subject_lonely_patch_indices) ==1
                        else
                            templatepatchtosubjectpatch_distance_mat{j,k}=distances(template_lonely_patch_indices,subject_lonely_patch_indices);
                        end
                    end
                end
            end
        end
        if exist('templatepatchtosubjectpatch_distance_mat','var')==0
        else
            templatepatchtosubjectpatch_distance_mat_by_net{i,1}=templatepatchtosubjectpatch_distance_mat;
            clear templatepatchtosubjectpatch_distance_mat
        end
    end
    save([output_subject_path filesep output_file_name '_templatepatchto' output_file_name 'patch_distance_mat_by_net.mat'],'templatepatchtosubjectpatch_distance_mat_by_net') %save this for later just in case.
end



for i=1:size(templatepatchtosubjectpatch_distance_mat_by_net,1)
    if isempty(templatepatchtosubjectpatch_distance_mat_by_net{i,1}) == 1
    else
        min_distance_of_every_grey_per_patch_size_checked = double([]);
        min_distance_of_every_grey_per_patch = double([]);
        for j=1:size(templatepatchtosubjectpatch_distance_mat_by_net{i,1},2)
            for k=1:size(templatepatchtosubjectpatch_distance_mat_by_net{i,1},1)
                min_distance_of_every_grey_per_patch(k,j) = min(min(templatepatchtosubjectpatch_distance_mat_by_net{i,1}{k,j}));
                %Do cluster check based on size also.
                if size(templatepatchtosubjectpatch_distance_mat_by_net{i,1}{k,j},1) < min_num_of_grays || size(templatepatchtosubjectpatch_distance_mat_by_net{i,1}{k,j},2) < min_num_of_grays
                    min_distance_of_every_grey_per_patch_size_checked(k,j) = NaN;
                else
                    min_distance_of_every_grey_per_patch_size_checked(k,j) = min(min(templatepatchtosubjectpatch_distance_mat_by_net{i,1}{k,j}));
                end
            end
        end
        min_distance_of_every_grey_per_patch_size_checked_allnets{i,1} = min_distance_of_every_grey_per_patch_size_checked;
        min_distance_of_every_grey_per_patch_allnets{i,1} = min_distance_of_every_grey_per_patch;
        clear min_distance_of_every_grey_per_patch
        clear min_distance_of_every_grey_per_patch_size_checked
    end
end

%Now that we've calculated the minimum distances, let's assign the
%networks, shall we?

%make copies of these for later. Because they will be modified to
%remove patches from the pool. (Akin to sample without replacement)
template_unmatched_clusters_poststep1_mod = template_unmatched_clusters_poststep1;
assignment_mat_post2 =assignment_mat; %save for later;

for n = 1: size(min_distance_of_every_grey_per_patch_allnets,1)
    if n == 4 || n ==6
    else
        min_distance_of_every_grey_per_patch = min_distance_of_every_grey_per_patch_allnets{n,1};
        possible_patch_values_template2 = 1:size(template_unmatched_clusters_poststep1_mod{n,1},1);
        
        for i = 1:size(min_distance_of_every_grey_per_patch,1)
            minimum = min(min(min_distance_of_every_grey_per_patch(min_distance_of_every_grey_per_patch<30))); % this is yucky, but essentially "Find the minimum value that is less than 30".
            if isempty(minimum) ==1
                disp('No more unmatched template patches found within 30 mm of an unmatched subject patch.')
                break
            else
                [x,y]=find(min_distance_of_every_grey_per_patch ==minimum);
                if length(x) > 1 %if there are multiple minimum values,
                    [minx,yidx] = min(x);
                    x=minx;
                    y=y(yidx);
                end
                min_distance_of_every_grey_per_patch(x,:)=nan;
                min_distance_of_every_grey_per_patch(:,y)=nan;
                
                %assign nets
                assignment_mat_post2{n,1}(subject_unassigned_patches{n,1}(y,1))= template_unmatched_clusters_poststep1_mod{n,1}(i);
                template_patches_to_exclude_post1(i) = x;
                subject_patches_no_longer_excluded_post1(i) =y;
            end
        end
        if exist('template_patches_to_exclude_post1', 'var') ==1
            subject_unassigned_patches{n,1}(subject_patches_no_longer_excluded_post1)=[];
            template_unmatched_clusters_poststep1_mod{n,1}(template_patches_to_exclude_post1) = [];
            clear subject_patches_no_longer_excluded_post1 template_patches_to_exclude_post1
        end
    end
end


%% Next steps: probably make dscalars to check that step 2 worked the way it
%should, and provide some stats.

for n=1:size(assignment_mat_post2,1)
    unlabeled_patches_post2(n,1) = size(assignment_mat_post2{n,1},1)-nnz(assignment_mat_post2{n,1});
    labeled_patches_post2(n,1)= nnz(assignment_mat_post2{n,1});
end

subject_patch_matched_matrix_post2=zeros(size(subject_patch_matrix,1),size(subject_patch_matrix,2));

for i= 1:size(assignment_mat_post2,1)
    for j=1:size(assignment_mat_post2{i,1},1)
        orig_patch_dix=subject_patch_matrix(:,i)==j;
        subject_patch_matched_matrix_post2(orig_patch_dix,i)=assignment_mat_post2{i,1}(j);
    end
end

if save_matched_dscalars==1
    %k=1;
    for net_num =net_list
        disp('Saving distance-checked patch-matched dscalars for each network for subject...')
        %subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '.dscalar.nii'];
        patch_matched_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched_dstance_' num2str(min_dist) '.dscalar.nii'];
        %subject_patch_matched_list{k,1} = patch_matched_subject_outputname_cifti_file;
        
        %if exist(patch_matched_subject_outputname_cifti_file,'file') ~=0
        %    disp('Subject patch dscalars has already been made for this network.')
        %else
        this_net_double = double(subject_patch_matched_matrix_post2(:,net_num));
        subject_cifti_obj.cdata=this_net_double;
        ciftisave(subject_cifti_obj,patch_matched_subject_outputname_cifti_file,wb_command);
        %system(cmd);
        %end
        %k=k+1;
    end
else
end

total_unlabelled_patches_post2 = sum(unlabeled_patches_post2);
total_labelled_patches_post2 = sum(labeled_patches_post2);
patch_assignedgreys_post2 = any(subject_patch_matched_matrix_post2,2);
number_of_grey_assigned_with_jaccard_only_post2 = sum(patch_assignedgreys_post2);
disp(['The number of subject unlabelled patched with jaccard only is: ' num2str(total_unlabelled_patches_post2) '/' num2str(max(subject_full_patches_list)) ]);
disp(['The number of  subject labelled patched with jaccard only is: ' num2str(total_labelled_patches_post2) '/' num2str(max(full_patches_list)) ' template patches ']);
disp(['The number of greys assigned with jaccard only is: ' num2str(number_of_grey_assigned_with_jaccard_only_post2) '/' num2str(size(patch_assignedgreys,1))]);
pause(4) % so that people can read
toc

disp('Continuing...')
%Next steps

disp('Next step: Assigning unassigned patches based on...')
disp('Step 3/3 - Distance to assigned clusters. Will have same assingment number as previosuly asigned.')

%% Start Step 3

i=1;k=1;
for i=1:size(assignment_mat_post2,1)
    %if isempty(assignment_mat{i,1}) ~=1
    subject_unassigned_patches_post2{i,1} = find(assignment_mat_post2{i,1}==0);
    %adjusted_missing_nets{i,1}=subject_list_of_adjusted_patches_by_networks{i,1}(unassigned_patches);
    for k=1:size(subject_unassigned_patches_post2{i,1},1)
        subject_adjusted_missing_nets_indices_by_net_post2{i,k}=find(subject_patch_matrix(:,i)==subject_unassigned_patches{i,1}(k));
    end
    num_clusters_unmatched_subject_clusters_by_net(i,1) =  length(assignment_mat_post2{i,1});
    %end
end

subject_unmatched_patch_indices_post2 = cell(size(assignment_mat_post2,1),max(num_clusters_unmatched_subject_clusters_by_net));
subject_matched_patch_indices_post2 = cell(size(assignment_mat_post2,1),max(num_clusters_unmatched_subject_clusters_by_net));
%This is redudant with the prior for loop but probably better since it
% get the indicies of all patches.
for i=1:size(assignment_mat_post2,1)
    if isempty(assignment_mat_post2{i,1}) ~=1
        for j=1:size(assignment_mat_post2{i,1},1)
            if assignment_mat_post2{i,1}(j) ==0
                subject_unmatched_patch_indices_post2{i,j} = find(subject_patch_matrix(:,i)==j);
            else
                subject_matched_patch_indices_post2{i,j} = find(subject_patch_matched_matrix_post2(:,i)==assignment_mat_post2{i,1}(j));
            end
        end
    end
end

%Now that we have the indicies for the missing and matched subject patches.
%get the distances.

tic
if exist('distances','var') ==1
    disp('Distance matrix already loaded.')
else
    disp('loading distance matrix...')
    if keep_cortical_subcortical_seperation ==1
        disp('Note: Distance matrix has eucliean distances between the cortex and subcortex set to 255mm (max uint8).')
        %load([support_folder filesep 'EUGEODistancematrix_XYZ_255interhem_unit8.mat'],'distances');
        load(distance_matrix_to_use,'distances');
    else
        disp('Note: Distance matrix uses eucliean distances between the cortex and subcortex.')
        %load([support_folder filesep 'EUGEODistancematrix_XYZ_unit8.mat'],'distances');
        load(distance_matrix_to_use,'distances');
    end
    toc
    
end

subjectunmatched2matched_distance_mat_by_net = cell(size(assignment_mat_post2,1),1);
for i=1:size(subject_unmatched_patch_indices_post2,1)
    %         if sum(~cellfun(@isempty,subjectunmatched2matched_distance_mat_by_net{i,:})>0) %check to see if there is naything to match (i.e. networks 4 || 6).
    %             subjectunmatched2matched_distance_mat_by_net{i,1}=[];
    %         else
    %         if i ==1 && isempty(subject_unmatched_patch_indices_post2{i,j}) ==1
    %             subjectpatchtosubjectpatch_distance_mat_by_net_post2{i,1}=[];
    %         else
    %num_unpatched_subject_patches=sum(~cellfun(@isempty,subject_unmatched_patch_indices_post2(i,:)));
    for j=1:size(subject_unmatched_patch_indices_post2,2)
        subject_lonely_patch_indices = subject_unmatched_patch_indices_post2{i,j};
        
        if isempty(subject_lonely_patch_indices) ==1
            subjectunmatched2matched_distance_mat{1,j} = [];
        else
            for k=1:size(subject_matched_patch_indices_post2,2)
                subject_matched_patch_indices = subject_matched_patch_indices_post2{i,k};
                if isempty(subject_matched_patch_indices) ==1
                    subjectunmatched2matched_distance_mat{k,j}=[];
                else
                    subjectunmatched2matched_distance_mat{k,j}=distances(subject_lonely_patch_indices,subject_matched_patch_indices);
                end
            end
        end
    end
    %end
    if exist('subjectunmatched2matched_distance_mat','var')==0
    else
        subjectunmatched2matched_distance_mat_by_net{i,1}=subjectunmatched2matched_distance_mat;
        clear subjectunmatched2matched_distance_mat
    end
    %end
end

%check to see if there is naything to match (i.e. networks 4 || 6).
for i=1:size(subjectunmatched2matched_distance_mat_by_net,1)
    if sum(sum(~cellfun(@isempty,subjectunmatched2matched_distance_mat_by_net{i,1})))>0
    else
        subjectunmatched2matched_distance_mat_by_net{i,1}=[];
    end
end

for i=1:size(subjectunmatched2matched_distance_mat_by_net,1)
    if isempty(subjectunmatched2matched_distance_mat_by_net{i,1}) == 1
    else
        %        min_distance_of_every_grey_per_patch_size_checked = double([]);
        %min_distance_of_every_grey_per_patch = cell(size(subjectunmatched2matched_distance_mat_by_net{i,1},1),size(subjectunmatched2matched_distance_mat_by_net{i,1},2));
        min_distance_of_every_grey_per_patch_nans=nan(size(subjectunmatched2matched_distance_mat_by_net{i,1},1),size(subjectunmatched2matched_distance_mat_by_net{i,1},2));
        for j=1:size(subjectunmatched2matched_distance_mat_by_net{i,1},2)
            for k=1:size(subjectunmatched2matched_distance_mat_by_net{i,1},1)
                if isempty(min(min(subjectunmatched2matched_distance_mat_by_net{i,1}{k,j}))) ==1
                    %there is nothing in the cell
                else
                    min_distance_of_every_grey_per_patch_nans(k,j) = min(min(subjectunmatched2matched_distance_mat_by_net{i,1}{k,j}));
                end
                %Do cluster check based on size also.
                %                 if size(subject_unmatched_patch_indices_post2{i,1}{k,j},1) < min_num_of_grays || size(subject_unmatched_patch_indices_post2{i,1}{k,j},2) < min_num_of_grays
                %                     min_distance_of_every_grey_per_patch_size_checked(k,j) = NaN;
                %                 else
                %                     min_distance_of_every_grey_per_patch_size_checked(k,j) = min(min(subject_unmatched_patch_indices_post2{i,1}{k,j}));
                %                 end
            end
        end
        %        min_distance_of_every_grey_per_patch_size_checked_allnets{i,1} = min_distance_of_every_grey_per_patch_size_checked;
        subject2subject_min_distance_of_every_grey_per_patch_allnets{i,1} = min_distance_of_every_grey_per_patch_nans;
        clear min_distance_of_every_grey_per_patch_nans
    end
end

assignment_mat_post3=assignment_mat_post2; % save for later.
subject_patch_matched_matrix_post3= subject_patch_matched_matrix_post2;
for n = 1: size(subject2subject_min_distance_of_every_grey_per_patch_allnets,1)
    if n == 4 || n ==6
    else
        min_distance_of_every_grey_per_patch_nans = subject2subject_min_distance_of_every_grey_per_patch_allnets{n,1};
        for i = 1:size(min_distance_of_every_grey_per_patch_nans,1)
            minimum = min(min(min_distance_of_every_grey_per_patch_nans)); % this is yucky, but essentially "Find the minimum value that is less than 30".
            if isempty(minimum) ==1 || isnan(minimum)
                disp('No more unmatched subject patches to match to known patches for this network.')
                break
            else
                [x,y]=find(min_distance_of_every_grey_per_patch_nans ==minimum); %y is old, x is new.
                if length(x) > 1 %if there are multiple minimum values,
                    [minx,yidx] = min(x);
                    x=minx;
                    y=y(yidx);
                end
                
                
                %assign nets
                %assignment_mat_post3{n,1}(subject_unassigned_patches_post2{n,1}(y,1))= template_matched_clusters_poststep1_mod{n,1}(i);
                %Be careful with the assignment number in the next line.
                %The value of x is index of the patch, not the "patch number."
                assignment_mat_post3{n,1}(y,1)= assignment_mat_post3{n,1}(x,1);
                
                subject_patch_matched_matrix_post3(subject_unmatched_patch_indices_post2{n,y},n) = assignment_mat_post3{n,1}(x,1);
                %min_distance_of_every_grey_per_patch(x,:)=nan;
                min_distance_of_every_grey_per_patch_nans(:,y)=nan;
                
                %subject_unmatched_patches_to_exclude_post1(i) = x;
                % do not remove patch numbers from the pool.
                %subject_patches_no_longer_excluded_post1(i) = y;
            end
        end
    end
end

%% Next steps: probably make dscalars to check that step 2 worked the way it
%should, and provide some stats.

for n=1:size(assignment_mat_post3,1)
    unlabeled_patches_post3(n,1) = size(assignment_mat_post3{n,1},1)-nnz(assignment_mat_post3{n,1});
    labeled_patches_post3(n,1)= nnz(assignment_mat_post3{n,1});
end

if save_matched_dscalars==1
    %k=1;
    for net_num =net_list
        disp('Saving distance-checked, patch-matched dscalars, proximity-matched for each network for subject...')
        %subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '.dscalar.nii'];
        step3_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched_dstance_' num2str(min_dist) '_post3.dscalar.nii'];
        %subject_patch_matched_list{k,1} = patch_matched_subject_outputname_cifti_file;
        
        if exist(step3_subject_outputname_cifti_file,'file') ~=0
            disp(['Subject patch dscalars has already been made for this network:' step3_subject_outputname_cifti_file])
        else
            this_net_double = double(subject_patch_matched_matrix_post3(:,net_num));
            subject_cifti_obj.cdata=this_net_double;
            ciftisave(subject_cifti_obj,step3_subject_outputname_cifti_file,wb_command);
            %system(cmd);
        end
        %k=k+1;
    end
else
end

total_unlabelled_patches_post3 = sum(unlabeled_patches_post3);
total_labelled_patches_post3 = sum(labeled_patches_post3);
patch_assignedgreys_post3 = any(subject_patch_matched_matrix_post3,2);
number_of_grey_assigned_with_jaccard_only_post3 = sum(patch_assignedgreys_post3);
disp(['The number of subject unlabelled patched with jaccard only is: ' num2str(total_unlabelled_patches_post3) '/' num2str(max(subject_full_patches_list)) ]);
disp(['The number of  subject labelled patched with jaccard only is: ' num2str(total_labelled_patches_post3) '/' num2str(max(full_patches_list)) ' template patches ']);
disp(['The number of greys assigned with jaccard only is: ' num2str(number_of_grey_assigned_with_jaccard_only_post3) '/' num2str(size(patch_assignedgreys_post3,1))]);
%pause(4) % so that people can read
toc

for i=1:size(patch_matrix,2)
    template_patches_this_net{i,1} = nonzeros(unique(patch_matrix(:,i)));
    subject_patches_this_net{i,1} = nonzeros(unique(subject_patch_matched_matrix_post3(:,i)));
end

subject_adjusted_patch_matrix_post3 = zeros(size(subject_patch_matched_matrix_post3,1),size(subject_patch_matched_matrix_post3,2));

for j = 1: size(subject_patch_matched_matrix_post3,2)
    if j ==1
        patches_this_net = nonzeros(unique(subject_patch_matched_matrix_post3(:,j)));
        subject_full_patches_list_post3 = patches_this_net;
        unique_patches_so_far = nonzeros(unique(subject_full_patches_list_post3));
        subject_list_of_adjusted_patches_by_networks_post3{j,1} = patches_this_net;
        for k = 1:size(patches_this_net,1)
            this_patches_indices = find(subject_patch_matched_matrix_post3(:,j)==patches_this_net(k));
            subject_adjusted_patch_matrix_post3(this_patches_indices,j) = patches_this_net(k);
        end
        
    else
        patches_this_net = nonzeros(unique(subject_patch_matched_matrix_post3(:,j)));
        if isempty(unique_patches_so_far) ==1
        else
            for k = 1:size(patches_this_net,1)
                this_patches_indices = find(subject_patch_matched_matrix_post3(:,j)==patches_this_net(k));
                subject_adjusted_patch_matrix_post3(this_patches_indices,j) = patches_this_net(k)+template_num_patches(j-1);
            end
            adjusted_patches_this_net = patches_this_net+template_num_patches(j-1);
            subject_list_of_adjusted_patches_by_networks_post3{j,1} = adjusted_patches_this_net;
            subject_full_patches_list_post3 = [subject_full_patches_list_post3; adjusted_patches_this_net];
        end
    end
end

subject_patch_label_vector =sum(subject_adjusted_patch_matrix_post3,2);
subject_unique_adjusted_patches = nonzeros(unique(subject_patch_label_vector));
template_unique_adjusted_patches = nonzeros(unique(template_patch_label_vector));

[unmatchable_log] = ~ismember(template_unique_adjusted_patches,subject_unique_adjusted_patches);
[matchable_log,matched_nets] = ismember(template_unique_adjusted_patches,subject_unique_adjusted_patches);
unmatchable_num = sum(unmatchable_log);
unmatched_indices = template_unique_adjusted_patches(unmatchable_log);
disp(['The number of patches that could not be matched is ' num2str(unmatchable_num) '/' num2str(max(template_unique_adjusted_patches)) ]);
disp(['The subject did not have the following networks: ' num2str(unmatched_indices')])
subject_cifti_obj.cdata=subject_patch_label_vector;
final_patch_path=[output_subject_path filesep output_file_name '_maxcombo' num2str(maximum_combination_of_nets) 'subject_net_all_unique_patches_post3.dscalar.nii'];
ciftisave(subject_cifti_obj,final_patch_path,wb_command);

disp('removing files')
for net_num =net_list
    subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_raw.dscalar.nii'];
    patch_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) 'patches_size' num2str(min_patch_size) 'patches.dscalar.nii'];
    step3_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched_dstance_' num2str(min_dist) '_post3.dscalar.nii'];
    patch_matched_subject_outputname_cifti_file = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched.dscalar.nii'];
    patch_matched_subject_outputname_cifti_file2 = [output_subject_path filesep output_file_name 'subject_net_' num2str(net_num) '_maxcombo' num2str(maximum_combination_of_nets) '_patch_matched_dstance_' num2str(min_dist) '.dscalar.nii'];
    cmd = sprintf('rm -v "%s"', subject_outputname_cifti_file);
    cmd2 = sprintf('rm -v "%s"', patch_subject_outputname_cifti_file);
    cmd3 = sprintf('rm -v "%s"', step3_subject_outputname_cifti_file);
    cmd4 = sprintf('rm -v "%s"', patch_matched_subject_outputname_cifti_file);
    cmd5 = sprintf('rm -v "%s"', patch_matched_subject_outputname_cifti_file2);
    system(cmd);
    system(cmd2);
    system(cmd3);
    system(cmd4);
    system(cmd5);
end
    
disp('done')


end