function [net_variance_mat, net_mean_mat, counts, mybins] = dconn_variance_per_network(dconn_file, assignments_vector_file,output_name,parcel_file)

%This function gets the variance of the dconn, after is's been sorted by
%each network, then gets the variance per network.

%R.Hermosillo 04/05/2021
%R.Hermosillo updated 02/25/2024 - updated to include make me option for parcel file.

%inputs are:
%dconn_files=path to a dense connectivity matrix cifti file (.dconni.nii).  This can also be a numeric matrix.

%assignments_vector_file=(string) path to a vector file that is the same length as the dconn file (.dscalar.nii).  This can also be a .csv. This can also be a numeric variable.

%outputname=(string) Some output file name. 

%parcel_file=path to a parcel file that contains the network assignments.  If one has
%not been created yet, and corresponds with a 14 network template matching
%solution, then you can write 'makeme' and the parcel will be automatically
%generated.


this_code = which('simple_cifti_average');
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
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/internal/utilities/MergeTimeSeries'));

warning('on')
wb_command=settings.path_wb_c; %path to wb_command

if isnumeric(dconn_file) ==1 
    dconn = dconn_file;
    clear dconn_file;
elseif strcmp(dconn_file(end-3:end),'.nii') ==1
    cii = ciftiopen(dconn_file,wb_command);
    dconn = cii.cdata;
    clear cii;
else
    load(dconn_file,'avg_cifti');
    dconn = avg_cifti;
    clear avg_cifti
end

%cii = ciftiopen('/panfs/roc/groups/3/rando149/shared/projects/ADHD_MedChal/template_matching_all_minutes/sub-1016_ses-both_merged_tasks_template_matched_Zscored_recolored.dscalar.nii',wb_command);
%assings = cii.cdata;
if isnumeric(assignments_vector_file)==1
    assigns_only_assigned = assignments_vector_file;
elseif strcmp(assignments_vector_file(end-3:end),'.nii') ==1
    assignscii=ciftiopen(assignments_vector_file,wb_command);
    assigns_only_assigned = assignscii.cdata;
elseif strcmp(assignments_vector_file(end-3:end),'.csv') ==1
    assigns_only_assigned_struct = importdata(assignments_vector_file);
    assigns_only_assigned =assigns_only_assigned_struct.data;
else
    load(assignments_vector_file,'assigns_only_assigned');
end

if exist([output_name '_TM_mean_per_network.mat'],'file') == 2
    disp(['Network .mat file found. It will not be remade: ' output_name '_TM_mean_per_network.mat'])
    load([output_name '_TM_mean_per_network.mat'], 'net_variance_mat', 'net_mean_mat','counts','mybins','countsnorm', 'net_variance_mat_alpha_sorted', 'net_mean_mat_alpha_sorted');
else
    
    
    %assigns = assigns_only_assigned;
     assigns = assigns_only_assigned(assigns_only_assigned>0);
   
    %assigns = assignments_vector;
    [sorted_networks,I] = sort(assigns); % get sorted indices;
    networks = unique(assigns); % get network assingments from template.
    disp('Sorting dconn...')
    sorted_dconn1 = dconn(I,I);
    %clear dconn newdconn
    
    net_start_idx = zeros(1, size(networks,1)); % preallocate for speed
    net_end_idx = zeros(1, size(networks,1)); % preallocate for speed
%     netRGBs = [
%         255 0 0;
%         0 0 153
%         255 255 0
%         0 255 0
%         13 133 160
%         50 50 50
%         102 0 204
%         102 255 255
%         255 128 0
%         178 102 153
%         0 102 153
%         102 255 102
%         60 60 251
%         200 200 200]/255;
    %matlab variance
    for i=1:size(networks,1)
        netidices = find(sorted_networks == networks(i));
        net_start_idx(i) = min(netidices)-1;
        net_end_idx (i)= max(netidices);
    end
    
    mybins = linspace(0, 1,501);
    net_variance_mat = zeros(size(networks,1),size(networks,1));
    net_mean_mat = zeros(size(networks,1),size(networks,1));
    counts = cell(size(networks,1),size(networks,1));
    countsnorm = cell(size(networks,1),size(networks,1));
    for i=1:size(networks,1)
        for j=1:size(networks,1)
            small_mat = sorted_dconn1(net_start_idx(i)+1:net_end_idx(i),net_start_idx(j)+1:net_end_idx(j));
            
            if i ==j %to get variance of within network connectivity, use upper triangle of matrix.
                small_matt = small_mat.';
                m_triulog = triu(true(size((small_matt))),1);
                net_net_mat_vec = small_matt(m_triulog).';
                
            else
                net_net_mat_vec = small_mat(:);
            end
            
            net_variance_mat(i,j) = var(net_net_mat_vec);
            net_mean_mat(i,j) = mean(net_net_mat_vec);
            h = histogram(net_net_mat_vec,mybins,'FaceColor','none','EdgeColor','r','DisplayStyle','stairs','LineWidth',2);
            counts{i,j} = h.Values;
            hnorm = histogram(net_net_mat_vec,mybins,'FaceColor','none','EdgeColor','k','DisplayStyle','stairs','LineWidth',2,'Normalization','probability');
            countsnorm{i,j} = hnorm.Values;
            %plot(mybins(1:end-1),counts{1,1})
            clear h net_net_mat_vec small_mat
            close all
        end
    end
    if strcmp(parcel_file,'makeme') ==1
        parcel = build_high_density_parcel_file(assigns,output_name);
    else
        load(parcel_file, 'parcel')
    end

    if isfield(parcel,'power_val')
        [~, alpha_sort]=sort([parcel.power_val],'ascend');
    else
        disp('networks are already sorted alphanumerically hopefully.')
        alpha_sort = 1:size([parcel.n],2);
    end
    disp('Saving mean and variance matrices sorted alphabetically too so they can match showM plots. ')
    net_variance_mat_alpha_sorted = net_variance_mat(alpha_sort,alpha_sort);
    net_mean_mat_alpha_sorted = net_mean_mat(alpha_sort,alpha_sort);
    %save([output_name '.mat'], 'net_variance_mat', 'net_mean_mat','counts','mybins')
    disp('Saving connectivity data:...' )
    save([output_name '_TM_mean_per_network.mat'], 'net_variance_mat', 'net_mean_mat','counts','mybins','countsnorm', 'net_variance_mat_alpha_sorted', 'net_mean_mat_alpha_sorted');
end

end