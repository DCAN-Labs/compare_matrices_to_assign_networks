function twins_mapping_wrapper(dt_or_ptseries_conc_file,motion_file,left_surface_file, right_surface_file, cifti_output_folder)
%BLV  maxNumCompThreads=2
%R. Hermosillo 1/8/2019
% this code takes in dtseries, motion, surfaces, for subject pairs and
% caluclates mtual information between individualized network assignments.
%
%BLV dt_or_ptseries_conc_file:

%hardcodes:
num_sub=26;
FD_threshold = 0.2;
smoothing_kernal = 1.75;
bit8 = 0;
TR=0.8;
minutes_limit = 'none';
series = 'dtseries';
data_type = 'dense';
wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 OMP_NUM_THREADS=2 /usr/local/bin/wb_command';
transform_data = 'Convert_FisherZ_to_r';
template_path = '/mnt/max/shared/code/internal/analyses/compare_matrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_fromsurfonly.mat';
remove_dconn =1;
make_dconn_conc = 0;

%% Start
%import concs
dtseries_file = importdata(dt_or_ptseries_conc_file);
motion_all = importdata(motion_file);
all_L_surface = importdata(left_surface_file);
all_R_surface = importdata(right_surface_file);

if remove_dconn ==1
    disp('Settings are set to remove all dconns after template matching.')
else
    disp('Settings are set to save all dconns. Be mindful of space.')
end


for i = 1:num_sub %number of subjects
    
    %step 0: infer output cifti name based on dtseries\
    conc = dtseries_file{i};
    motion = motion_all{i};
    L_surface = all_L_surface{i};
    R_surface = all_R_surface{i};
    
    [~,filename_long] = fileparts(conc);
    [~,filename_short]= fileparts(filename_long);
    output_cifti_name =[filename_short '_template_matched'];
    
    %Step 1: make matrix
    %BLV added additional path for template_matching_RH
    support_folder=['/mnt/max/shared/code/internal/analyses/compare_matrices/'];
    addpath(genpath(support_folder));
    addpath('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/');
    subjectdconn = cifti_conn_matrix(conc,series,motion, FD_threshold, TR, minutes_limit, smoothing_kernal,L_surface,R_surface,bit8,make_dconn_conc);
    
    %Step 2: get network assingments
    [~, eta_subject_index{i}] = template_matching_RH(subjectdconn, data_type, template_path,transform_data,output_cifti_name, cifti_output_folder ,wb_command);
    
    if remove_dconn ==1 % RH added in case filespace becomes an limited.
        % Step 2.5: Remove dconn to save space.
        cmd = ['rm -f ' subjectdconn];
        disp(cmd);
        system(cmd);
    else
        disp('Keeping dconn. Be mindful of space.')
    end
end

disp('Done calculating network assingments for all pairs. Proceeding to calculate mutual information for all pairs.');



for i = 1:2:num_sub %number of subjects
    %Step 3: calculate mutualinformation
    muI(i,1) = MutualInformation(eta_subject_index{i},eta_subject_index{i+1}); %Mutual information
    [VIn(i,1), MIn(i,1)] = partition_distance(eta_subject_index{i},eta_subject_index{i+1}); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
end

%display summary
MIn
MuI
VIn

dlmwrite([cifti_output_folder '/' output_cifti_name '.txt'],MIn,muI,VIn,'roffset',1)

disp('done calculating mutual information for all pairs')

end