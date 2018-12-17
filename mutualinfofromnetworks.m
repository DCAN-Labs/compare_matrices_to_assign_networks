function mutualinfofromnetworks(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_vector,include_all_frames, smoothing_kernal,left_surface_file, right_surface_file, bit8, output_cifti_name,method, cifti_enhancement,other_half_networks_cii)
%% This code is designed to make correlation matrix from a subject's dtseries, and motion, and surface files (if smoothing is desired) See documentation for cifti_conn_matrix.
%This code then runs template matching on subjects.  It adds

if isnumeric(minutes_vector)==1
else
    if strcmp(minutes_vector,'none') == 1      
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
    

support_folder=[pwd '/support_files'];
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

warning('off') %supress addpath warnings to nonfolders.
disp('Attempting to add neccesaary paths and functions.')
for i=1:np
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
wb_command=settings.path_wb_c; %path to wb_command
warning('on')

tic
for i=1:length(minutes_vector)
    if strcmp(minutes_vector{i},'none') == 1
        minutes_name = 'none';
    else
        minutes_name = num2str(minutes_vector{i});
    end
    
    switch series
        case 'ptseries'
            
            data_type = 'parcellated';
            if cifti_enhancement ==1
                output_cifti_scalar_name  = [output_cifti_name minutes_name '_method_' method '_NE.pscalar.nii'];
            else
                output_cifti_scalar_name  = [output_cifti_name minutes_name '_method_' method '.pscalar.nii'];
            end
        case 'dtseries'
            data_type = 'dense';
            
            if cifti_enhancement ==1
                output_cifti_scalar_name  = [output_cifti_name minutes_name '_method_' method '_NE.dscalar.nii'];
            else
                output_cifti_scalar_name  = [output_cifti_name minutes_name '_method_' method '.dscalar.nii'];
            end
            
    end
    if exist(output_cifti_scalar_name) ==0
        disp([output_cifti_scalar_name ' not found. Making matrix prior to template matching.']);
        %% Step 1: make connectivity matrix
        addpath('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien');
        temp_name = cifti_conn_matrix(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_vector{i}, smoothing_kernal,left_surface_file, right_surface_file, bit8);
        %% step 1b: refine matrix? reduce noise in matrix.
        if cifti_enhancement ==1
        enhanced_network_name = run_cifti_network(temp_name,data_type);
        temp_name = enhanced_network_name;
        else
        end
        
        %% Step 2: Do template matching to get network assingments
        [new_subject_labels, output_cifti_scalar_name] = comparematrices_test(temp_name,[output_cifti_name minutes_vector{i}],method,data_type,cifti_enhancement);
        %% Step 3: Remove dconn.
        %conn_dir = fileparts(dt_or_ptseries_conc_file);
        cmd = ['rm -f ' temp_name];
        disp(cmd);
    else
        disp([output_cifti_scalar_name ' found. Loading scalar.']);
        current_cifti_labels = ciftiopen(output_cifti_scalar_name,wb_command);
        new_subject_labels = current_cifti_labels.cdata;
    end
    %% Step 4: load template networks from 2nd half of data:
    otherhalf_cii = ciftiopen(other_half_networks_cii,wb_command);
    other_half_networks = otherhalf_cii.cdata;
    
    %% Step 5: Calculate mutual information
    muI(i,1) = MutualInformation(other_half_networks,new_subject_labels); %Mutual information
    [VIn(i,1), MIn(i,1)] = partition_distance(other_half_networks,new_subject_labels); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
%
    correl(i,1) = corr(other_half_networks,new_subject_labels);
    toc
end





save([output_cifti_name 'MuI.mat'],'muI','correl','VIn','MIn');
disp(['Done calulating mutual information for all ' num2str(length(minutes_vector)) ' matrices.']);

end