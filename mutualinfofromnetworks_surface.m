function [muI] = mutualinfofromnetworks_surface(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR,include_all_frames, smoothing_kernal,left_surface_file, right_surface_file, bit8, output_cifti_name,community_detection, method, cifti_enhancement,other_half_networks_cii,num_interval_reps,indepen_time_series, additional_mask)
%% This code then runs template matching or infomap on subjects and calculate the mututal information to a specifided dscalar.
%%This code is designed to make correlation matrix from a subject's dtseries, and motion, and surface files (if smoothing is desired) See documentation for cifti_conn_matrix.
% Correlation matrices are automatically Z-scored.

%Hermosillo R. 1/25/2019
%--------------

%UPDATED 2/12/2019
%Updated to included option to run multiple repetitions 
%Update also one to sample the same dconn many times (in case you're interested in checking the stochatstic nature of infomap).

%How to run code
% Arguments are
% 1) dtseries = desnse time series
% 2) series = 'ptseries' or dtseries
% 3) motion_file = Power 2014 motion.mat file
% 4) FD_threshold = your frame-wise displament threshold. (e.g. the frames you want to keep from your time series.
% 5) TR = repitition time of your data.
% 6) minutes_vector = a cell vector that contain each desired minute.  Leave empty to specifiy 'none minutes limit' (i.e. make an 'all frames' dconn)
% 7) include_all_frames, set to 1 to make a include minutes limit in your minutes vector
% 8) smoothing_kernal = specify smoothing kernal (note: if no smoothing file to be used, type 'none')
% 9) left_surface_file = full path to the left midthickness surface file.
% 10) right_surface_file full path to the right midthickness surface file.
% 11) bit8 = set to 1 if you want to make the outputs smaller in 8bit. Not recommended.  Set this to 0.
% 12) output_cifti_name = The name of your output cifti.  This option only works for template matching.  It does not work for infomap.
% 13) community_detection =  provide either 'template_matching' or 'Infomap'
% 14) method = Only used in template matching.  Specifiy 'Template matching'.
% 15) cifti_enhancement  = Not debugged.  Set to 0.
% 16) other_half_networks_cii,
% 17) num_interval_reps.  Number of repetitions of community detection. This is different that the specified number of reps that info map runs. Each repetition will make it's own dscalar.nii
% 18) indepen_time_series = if 1, each repititon will use a different dconn, if 0 it will use the same dconn. When minutes are specified these are randomly sampled from available frames below the FD threshold to make the dconn.  If you want to use the exact same dconn, this will make symlinks to that dconn. 


%% Inititalize by adding paths
this_code = which('mutualinfofromnetworks_surface');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,10);

%load('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/minutes_vector.mat') %[1,2,3,4,5,10,15,20];

minutes_vector = [1,2,3,4,5,10,15,20];

if isnumeric(minutes_vector)==1
else
    if strcmp(minutes_vector,'none') == 1
        minutes_vector = cellstr('none');
    else
        minutes_vector = str2double(minutes_vector);
    end
end

if isempty(minutes_vector) %determine if minutes vector includes a 'none' option.
    if include_all_frames == 1
        minutes_vector = cellstr('none');
    else
        disp('No minutes vector and no "option to include all frames (1/0)" provided.  Check your input arguments.')
        return
    end
else
    minutes_vector = num2cell(minutes_vector);
    nonecell = cellstr('none');
    minutes_vector = [minutes_vector nonecell]; % concatenate cells with numbers and cell with 'none'
end

if isnumeric(num_interval_reps)==1
else
    num_interval_reps = str2double(num_interval_reps);
end

if isnumeric(indepen_time_series)==1
else
    indepen_time_series = str2double(indepen_time_series);
end

if isnumeric(cifti_enhancement)==1
else
    cifti_enhancement = str2double(cifti_enhancement);
end




warning('off') %supress addpath warnings to nonfolders.
disp('Attempting to add neccesaary paths and functions.')
for i=1:np
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase

addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/Matlab_CIFTI'))
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/gifti-1.6'))

wb_command=settings.path_wb_c; %path to wb_command
warning('on')
num_reps_vector = 1:num_interval_reps;


%some hardcodes for infomap:
template = 'none';
min_distance = 20;
tie_density = 0.02;
min_network_size = 400;
min_region_size = 30;
%community_detection = 'infomap';
num_reps = 20;
donotZscore = 0; 
surface_only = 1;

%hardcodes for template matching;
allow_overlap = 1; 
overlap_method = 'smooth_then_derivative';

%% Start

if indepen_time_series == 0 %sample from the same dconn.
   
    for r=1:num_interval_reps
        tic
        for i=1:length(minutes_vector)
            if strcmp(minutes_vector{i},'none') == 1
                minutes_andreps_name = 'none';
            else
                minutes_andreps_name = [num2str(minutes_vector{i}) 'reps' num2str(num_reps_vector(r))];
            end
            
                    data_type = 'dense';
                    %if cifti_enhancement ==1
                    %    output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '_NE.dscalar.nii'];
                    %else
                    templ_dscalar_out_dir = pwd;
                    if donotZscore == 0 
                    template_output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '_Zscored_recolored.dscalar.nii'];
                    else
                    template_output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '_recolored.dscalar.nii'];    
                    end
                    %end

                    
                    
            if  exist([templ_dscalar_out_dir '/' template_output_cifti_scalar_name],'file') == 0
                disp([templ_dscalar_out_dir '/' template_output_cifti_scalar_name ' not found. Making matrix prior to template matching.']);
                %% Step 1: make connectivity matrix
                %%HARDCODE WARNING
                addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/hcp_comm_det_damien');
                remove_outliers =0; %hardcoded options for making connectivity matrices. changed this to 0 for infants
                %additional_mask = 'none';
                orig_temp_name = cifti_conn_matrix(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_vector{i}, smoothing_kernal,left_surface_file, right_surface_file, bit8, remove_outliers,additional_mask);
                temp_short = orig_temp_name(1:end-10); %remove dconn.nii
                
                %Build symlinks first
                cmd = ['ln -sf ' orig_temp_name ' ' temp_short 'rep' num2str(r) '.dconn.nii'];
                system(cmd);
                disp(cmd)
                temp_name = [temp_short 'rep' num2str(r) '.dconn.nii'];
                Ztemp_name = [temp_short 'rep' num2str(r) 'Zscored.dconn.nii'];
                [~,file_Ztemp_name] = fileparts(Ztemp_name);
                
                 dir1 = [pwd '/Community_Detection_Min_Dist_' num2str(min_distance) '_TieDen_' num2str(tie_density) '_MinNet_Size_' num2str(min_network_size) '_MinReg_Size_' num2str(min_region_size)];
                 info_dscalar_out_dir = [dir1 '/' file_Ztemp_name '.nii/wb_ready_files'];
                 infomap_output_cifti_scalar_name  = [file_Ztemp_name(1:end-6) '_pajek_' num2str(tie_density*100) 'perc_' num2str(min_distance) 'dist_bin_No_COMMUNITIES_recolored.dscalar.nii'];

                %% step 1b: refine matrix? reduce noise in matrix.
%                 if cifti_enhancement ==1
%                     enhanced_network_name = run_cifti_network(temp_name,data_type);
%                     temp_name = enhanced_network_name;
%                 else
%                 end
                
                %% Step 2: Run community detection to get network assingments
                %switch community_detection
                        
                    %case 'infomap'
                        disp('Community detection method is infomap')
                        disp('Adding paths for infomap.')
                        addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair')
                        %HARDCODES
                        path_dist_matrix = '/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
                        
                        disp(['running community detection with command: cifti_community_detection_surface(' temp_name ',' path_dist_matrix ',' template ',' num2str(min_distance) ',' num2str(tie_density) ',' num2str(min_network_size) ',' num2str(min_region_size) ',' num2str(donotZscore) ',' num2str(surface_only) ',' 'infomap' ',' num2str(num_reps) ')']);
                        [ info_output_cifti_scalar_name] = cifti_community_detection_surface(temp_name,path_dist_matrix,template,min_distance,tie_density,min_network_size,min_region_size,donotZscore,surface_only,'infomap',num_reps);
                        cleaned_infomapfile = ciftiopen(info_output_cifti_scalar_name,wb_command);
                        new_subject_labels = cleaned_infomapfile.cdata;
                        
                     %case 'template_matching'
                        disp('Community detection method is template_matching')
                        [new_subject_labels, template_output_cifti_scalar_name] = comparematrices_test_surface(temp_name,[output_cifti_name minutes_andreps_name],method,data_type,cifti_enhancement,allow_overlap,overlap_method);
                       
                    %case 'bigclam'
%                         disp('big clam not yet supported')
%                         return
%                         template = 'none';
%                         min_distance = 20;
%                         tie_density = 0.015;
%                         min_network_size = 400;
%                         min_region_size = 30;
%                         %community_detection = 'Bigclam';
%                         num_reps = 20;
%                         
%                         disp(['Community detection method is Bigclam'])
%                         [new_subject_labels, output_cifti_scalar_name] = cifti_community_detection_surface(temp_name,distance_matrix,template,min_distance,tie_density,min_network_size,min_region_size,community_detection,num_reps);
%                     otherwise
                %end
                %% Step 3: Remove dconn.
                %conn_dir = fileparts(dt_or_ptseries_conc_file);
                if strcmp(minutes_vector{i},'none') == 1 % do not remove dconns for all frames dconns as it is unecessary.
                else
                    cmd = ['rm -f ' orig_temp_name];
                    disp(cmd);
                    system(cmd);
                    cmd = ['rm -f ' Ztemp_name];
                    disp(cmd);
                    system(cmd);                    
                    
                    
                end
            else
                disp([templ_dscalar_out_dir '/' template_output_cifti_scalar_name ' found. Loading scalar.']);
                current_cifti_labels = ciftiopen(template_output_cifti_scalar_name,wb_command);
                new_subject_labels = current_cifti_labels.cdata;
            end
            %% Step 4: load template networks from 2nd half of data:
            otherhalf_cii = ciftiopen(other_half_networks_cii,wb_command);
            other_half_networks = otherhalf_cii.cdata;

            %% Step 5: Calculate mutual information
            
            disp(['Test scalar =' templ_dscalar_out_dir '/' template_output_cifti_scalar_name])            
            disp(['Reference scalar = ' other_half_networks_cii])
            disp('Caluating mutual information between scalars')
            %% Step 5: Calculate mutual information
            muI(i,r,1) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
            [VIn(i,r,1), MIn(i,r,1)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            %correl(i,r,1) = corr(other_half_networks,new_subject_labels);
            toc
  
        end
    end
    
else %sample from the same dconn.
    for i=1:length(minutes_vector)
        %% Step 1: make connectivity matrix
        %%HARDCODE WARNING
        addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/hcp_comm_det_damien');
        remove_outliers =0; %hardcoded options for making connectivity matrices. changed this to 0 for infants
        orig_temp_name = cifti_conn_matrix(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_vector{i}, smoothing_kernal,left_surface_file, right_surface_file, bit8, remove_outliers,additional_mask);
        temp_short = orig_temp_name(1:end-10); %remove dconn.nii
        [~,temp_file] = fileparts(orig_temp_name);
        %temp_file_short = temp_file(1:end-6);
        for r = 1:num_interval_reps
            
            %Build symlinks first
            cmd = ['ln -sf ' orig_temp_name ' ' temp_short 'rep' num2str(r) '.dconn.nii'];
            system(cmd);
            disp(cmd)
            temp_name = [temp_short 'rep' num2str(r) '.dconn.nii'];
            Ztemp_name = [temp_short 'rep' num2str(r) 'Zscored.dconn.nii'];
            [~,file_Ztemp_name] = fileparts(Ztemp_name);

            %Build output names for template matching
            if strcmp(minutes_vector{i},'none') == 1
                minutes_andreps_name = 'none';
            else
                minutes_andreps_name = [num2str(minutes_vector{i}) 'reps' num2str(num_reps_vector(r))];
            end
            
            switch series
                case 'ptseries'
                    data_type = 'parcellated';
                    if cifti_enhancement ==1
                        output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '_NE.pscalar.nii'];
                    else
                        output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '.pscalar.nii'];
                    end
                case 'dtseries'
                    data_type = 'dense';
                    if cifti_enhancement ==1
                        output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '_NE.dscalar.nii'];
                    else
                        output_cifti_scalar_name  = [output_cifti_name minutes_andreps_name '_method_' method '.dscalar.nii'];
                    end
                otherwise
                    disp('series type must either "ptseries" or "dtseries". check your inputs.')
                    return
            end
            
            switch community_detection
                case 'template_matching'
                    dscalar_out_dir = pwd;
                case 'infomap'
                    dir1 = [pwd '/Community_Detection_Min_Dist_' num2str(min_distance) '_TieDen_' num2str(tie_density) '_MinNet_Size_' num2str(min_network_size) '_MinReg_Size_' num2str(min_region_size)];
                    dscalar_out_dir = [dir1 '/' file_Ztemp_name '.nii/wb_ready_files'];
                    output_cifti_scalar_name  = [file_Ztemp_name(1:end-6) '_pajek_' num2str(tie_density*100) 'perc_' num2str(min_distance) 'dist_bin_No_COMMUNITIES_recolored.dscalar.nii'];
                case 'bigclam'
                otherwise        
            end
            
            if exist([dscalar_out_dir '/' output_cifti_scalar_name],'file') == 0
                disp([dscalar_out_dir '/' output_cifti_scalar_name ' not found. Making matrix prior to template matching.']);
                %% step 1b: refine matrix? reduce noise in matrix.
                if cifti_enhancement ==1
                    enhanced_network_name = run_cifti_network(temp_name,data_type);
                    temp_name = enhanced_network_name;
                else
                end
                
                %% Step 2: Run community detection to get network assingments
                
                switch community_detection
                    case 'template_matching'
                        disp('Community detection method is template_matching')
                        [new_subject_labels, output_cifti_scalar_name] = comparematrices_test_surface(temp_name,[output_cifti_name num2str(minutes_vector{i})],method,data_type,cifti_enhancement);
                        
                    case 'infomap'
                        disp('Community detection method is infomap')
                        disp('Adding paths for infomap.')
                        addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair')
                        %HARDCODES
                        
                        path_dist_matrix = '/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat';
                        disp(['running community detection with command: cifti_community_detection_surface(' temp_name ',' path_dist_matrix ',' template ',' num2str(min_distance) ',' num2str(tie_density) ',' num2str(min_network_size) ',' num2str(min_region_size) ',' num2str(donotZscore) ',' num2str(surface_only) ',' 'infomap' ',' num2str(num_reps) ')']);
                        [ output_cifti_scalar_name] = cifti_community_detection_surface(temp_name,path_dist_matrix,template,min_distance,tie_density,min_network_size,min_region_size,donotZscore,surface_only,'infomap',num_reps);
                        cleaned_infomapfile = ciftiopen(output_cifti_scalar_name,wb_command);
                        new_subject_labels = cleaned_infomapfile.cdata;
                        
                    case 'bigclam'
                        disp('big clam not yet supported')
                        return
                        template = 'none';
                        min_distance = 20;
                        tie_density = 0.02;
                        min_network_size = 400;
                        min_region_size = 30;
                        %community_detection = 'Bigclam';
                        num_reps = 20;
                        
                        disp(['Community detection method is Bigclam'])
                        [new_subject_labels, output_cifti_scalar_name] = cifti_community_detection_surface(temp_name,path_dist_matrix,template,min_distance,tie_density,min_network_size,min_region_size,donotZscore,surface_only,community_detection,num_reps);
                        
                    otherwise
                end
                %% Step 3: Remove dconn.
                %conn_dir = fileparts(dt_or_ptseries_conc_file);
                if strcmp(minutes_vector{i},'none') == 1 % do not remove dconns for all frames dconns as it is unecessary.
                else
                    cmd = ['rm -f ' temp_name];
                    disp(cmd);
                    system(cmd);
                end
            else
                disp([dscalar_out_dir '/' output_cifti_scalar_name ' found. Loading scalar.']);
                current_cifti_labels = ciftiopen([dscalar_out_dir '/' output_cifti_scalar_name],wb_command);
                new_subject_labels = current_cifti_labels.cdata;
            end
            %% Step 4: load template networks from 2nd half of data:
            otherhalf_cii = ciftiopen(other_half_networks_cii,wb_command);
            other_half_networks = otherhalf_cii.cdata;
            disp(['Test scalar = ' dscalar_out_dir '/' output_cifti_scalar_name])
            disp(['Reference scalar = ' other_half_networks_cii])
            disp('Caluating mutual information between scalars')
            %% Step 5: Calculate mutual information
            muI(i,r,1) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
            [VIn(i,r,1), MIn(i,r,1)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            %correl(i,r,1) = corr(other_half_networks,new_subject_labels);
            toc
        end
    end
end

save([output_cifti_name 'MuI.mat'],'muI','VIn','MIn','-v7.3');
disp(['Done calulating mutual information for all ' num2str(length(minutes_vector)) ' matrices.']);


end