function twins_mapping_wrapper(dt_or_ptseries_conc_file,motion_file,left_surface_file, right_surface_file, output_file_name, cifti_output_folder,TR,minutes_limit,FD_threshold,transform_data,template_path,surface_only,already_surface_only,use_all_ABCD_tasks, run_infomap_too,output_directory, dtseries_conc,use_continous_minutes,memory_limit_value)
%R. Hermosillo 1/8/2019
% this code takes in dtseries, motion, surfaces, for subject pairs and
% caluclates mtual information between individualized network assignments.
%

%hardcodes:
num_sub=length(dt_or_ptseries_conc_file);
%FD_threshold = 0.2;
smoothing_kernal = 2.55;
bit8 = 0;
%TR=2.5;
%minutes_limit = 'none';
series = 'dtseries';
data_type = 'dense';
%wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 OMP_NUM_THREADS=2 /usr/local/bin/wb_command';
%wb_command = '/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-1.2.3-HCP/bin_rh_linux64/wb_command';
%transform_data = 'Convert_FisherZ_to_r';
%transform_data = 'Convert_to_Zscores';
%template_path = '/mnt/max/shared/code/internal/analyses/compare_matrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_fromsurfonly.mat';
%template_path = '/mnt/max/shared/code/internal/analyses/compare_matclearrices/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat';
%template_path = '/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/seedmaps_ABCD164template_dtseries_all_networksZscored.mat';
%template_path = '/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/seedmaps_ABCD164template_SMOOTHED_dtseries_all_networksZscored.mat';
clean_up_intermed_files = 1;
make_dconn_conc = 0;
%output_file_name = 'ADHD315';
calculate_mutual_info = 0;
make_cifti_from_results = 1;
allow_overlap = 1;
overlap_method = 'smooth_then_derivative';
remove_outliers= 1; additional_mask = 'none';

%check input format
if isnumeric(TR)==1
else
    TR=str2num(TR);
end

if isnumeric(FD_threshold)==1
else
    FD_threshold=str2num(FD_threshold);
end

if isnumeric(surface_only)==1
else
    surface_only=str2num(surface_only);
end

if isnumeric(already_surface_only)==1
else
    already_surface_only=str2num(already_surface_only);
end

if isnumeric(use_all_ABCD_tasks)==1
else
    use_all_ABCD_tasks=str2num(use_all_ABCD_tasks);
end

if isnumeric(run_infomap_too)==1
else
    run_infomap_too=str2num(run_infomap_too);
end

if isnumeric(use_continous_minutes) ==1
else
    use_continous_minutes = str2num(use_continous_minutes);
end

if isnumeric(memory_limit_value) ==1
else
    memory_limit_value = str2num(memory_limit_value);
end


%Check to make sure that  minutes limit is a number (unless you've set it
%to 'none')
if strcmp(minutes_limit,'none')==1 || strcmp(minutes_limit,'None')==1 || strcmp(minutes_limit,'NONE')==1
    minutes_limit= 'none';
elseif isnumeric(minutes_limit)==1
    disp('minutes limit is passed in as numeric.')
else
    disp('minutes limit is passed in as string. converting to numeric')
    minutes_limit = str2num(minutes_limit);
end

%% Adding paths for this function
this_code = which('twins_mapping_wrapper');
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
rmpath('/home/faird/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/internal/utilities/MergeTimeSeries'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command


[dtpath, dtfile]=fileparts(dt_or_ptseries_conc_file);
split_dtpath=strsplit(dtfile,'_');
subID=[split_dtpath{1} '_' split_dtpath{2}];

if use_all_ABCD_tasks == 1
    % NOTE: this option has not been tested for conc files.
    [dt_conc_name,motion_conc_names] = make_scan_conc(dtpath,dtfile); %use dtseries file name and location to find other tasks.
    MergeTimeSeries('TimeSeriesFiles',dt_conc_name,'MotionFiles',motion_conc_names,'OutputFile',[dtpath filesep subID '_merged_tasks.dtseries.nii'],'MotionOutputFile',[dtpath filesep subID '_merged_tasks_motion.mat'])
    %Set time series and motion files to the newly merged data.
    dt_or_ptseries_conc_file= [dtpath filesep subID '_merged_tasks.dtseries.nii'];
    motion_file=[dtpath filesep subID '_merged_tasks_motion.mat'];
else
end

%% Start
%import concs

conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dtseries_file = importdata(dt_or_ptseries_conc_file);
    motion_all = importdata(motion_file);
    all_L_surface = importdata(left_surface_file);
    all_R_surface = importdata(right_surface_file);
    
else
    dtseries_file = {dt_or_ptseries_conc_file};
    motion_all = {motion_file};
    all_L_surface = {left_surface_file};
    all_R_surface = {right_surface_file};
end

%dtseries_file = importdata(dt_or_ptseries_conc_file);
%motion_all = importdata(motion_file);
%all_L_surface = importdata(left_surface_file);
%all_R_surface = importdata(right_surface_file);


for i = 1:length(dtseries_file)
    if exist(dtseries_file{i},'file') == 0
        disp(['NOTE Subject Series ' num2str(i) ' does not exist'])
        return
    else
    end
end
disp('All series files exist continuing ...')

if exist(template_path,'file') == 0
    disp('NOTE template does not exist does not exist')
    return
else
end


if clean_up_intermed_files ==1
    disp('Settings are set to remove all dconns after template matching.')
else
    disp('Settings are set to save all dconns. Be mindful of space.')
end

for i = 1:length(dtseries_file) %number of subjects
    
    %step 0: infer output cifti name based on dtseries\
    subject_dt_series = dtseries_file{i};
    motion = motion_all{i};
    L_surface = all_L_surface{i};
    R_surface = all_R_surface{i};
    
    [path_to_orig_dtseries,filename_long] = fileparts(subject_dt_series);
    [~,filename_short]= fileparts(filename_long);
    output_cifti_name =[filename_short '_template_matched'];
    output_cifti_name_info =[filename_short '_infomap'];
    
    
    switch data_type
        case 'parcellated'
            output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.pscalar.nii' ];
        case 'dense'
            
            output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.dscalar.nii' ];
    end
    
    if strcmp(transform_data,'Convert_to_Zscores') == 1
        output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '_Zscored.dscalar.nii' ];
    else
    end
    %Step 1: make matrix
    %support_folder=['/mnt/max/shared/code/internal/analyses/compare_matrices/'];
    %addpath('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/');
    
    %support_folder='/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices'; -old
    support_folder='/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks';
    %addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/hcp_comm_det_damien/'); -old
    addpath('/home/faird/shared/code/internal/utilities/cifti_connectivity/src');
    addpath(support_folder);
    
    
    
    %     if exist([cifti_output_folder '/' output_cifti_name '.mat']) == 2
    %         disp('.mat file already created.  loading...');
    %         load([cifti_output_folder '/' output_cifti_name '.mat']);
    %         new_subject_labels = eta_subject_index;
    %         eta_net_assign{i} = eta_subject_index;
    %
    %         if exist([output_cifti_scalar_name],'file') == 0
    %             disp('saving file to cifti')
    %             saving_template =ciftiopen(settings.path{8}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
    %             saving_template.cdata = single(new_subject_labels);
    %             %addpath('/mnt/max/shared/code/internal/utilities/corr_pt_dt/support_files');
    %             disp('Saving new scalar')
    %             save(saving_template, [cifti_output_folder '/' output_cifti_name '.gii'],'ExternalFileBinary') %DF: This save not pointing to the right place fix!
    %             %save(saving_template, output_cifti_scalar_name, wb_command) %RH sifti save fix.
    %             disp('Converting scalar .gii to .nii')
    %             unix([wb_command ' -cifti-convert -from-gifti-ext ' cifti_output_folder '/' output_cifti_name '.gii ' output_cifti_scalar_name ]);
    %             disp('Removing .gii')
    %             unix(['rm -f ' cifti_output_folder '/' output_cifti_name '.gii']);
    %             unix(['rm -f ' cifti_output_folder '/' output_cifti_name '.dat']);
    %else
    %    disp(['Cifti file already exists. Data was loaded from corresponding .mat file.' ]);
    %end
    dotsloc = strfind(output_cifti_scalar_name,'.');
    sub_basename = output_cifti_scalar_name(1:(dotsloc(end-1)-1));
    recolored_name = [sub_basename '_recolored.dscalar.nii'];
    
    if run_infomap_too ==1
        if use_all_ABCD_tasks == 1
            expected_final_dtseries_file = [cifti_output_folder filesep subID '_merged_tasks_infomap_densities_recolored.dtseries.nii'];
        else
            expected_final_dtseries_file = [cifti_output_folder filesep subID '_infomap_densities_recolored.dtseries.nii'];
        end
            expected_info_dscalar_name = [expected_final_dtseries_file(1:end-13),'.dscalar.nii'];
    end
    
    
    %% final_outputs_found check 
    if run_infomap_too  ==1
        if exist(output_cifti_scalar_name,'file') == 2 && exist(expected_info_dscalar_name,'file') == 2
            disp(['infomap final dtseries found for this subject: ' expected_info_dscalar_name])
            disp(['Template_matching scalar found for this subject: ' output_cifti_scalar_name])
            if exist(recolored_name,'file') == 0
                disp('Cleaned template matching file not found. cleaning...')
                clean_dscalars_by_size(output_cifti_scalar_name,[],[],[],[],30,[],0,1)
            else
                disp(['No cleaning necessary: Cleaned networks file found: ' recolored_name])
            end
            final_outputs_found =1;
        else
            final_outputs_found =0;
        end
    else
        if exist(output_cifti_scalar_name,'file') == 2
            disp(['Template_matching scalar found for this subject: ' output_cifti_scalar_name])
            if exist(recolored_name,'file') == 0
                disp('Cleaned template matching file not found. cleaning...')
                clean_dscalars_by_size(output_cifti_scalar_name,[],[],[],[],30,[],0,1)
            else
                disp(['No cleaning necessary: Cleaned networks file found: ' recolored_name])
            end
            final_outputs_found =1;
        else
            final_outputs_found =0;
        end     
    end
    
    if final_outputs_found ==1
        disp('Final outputs found.  Lets skip a bunch of steps to the picture-making part.')
    else % start from beginning
        [motion_path, motion_name_only] =fileparts(motion);
        if surface_only ==1
            if  already_surface_only ==0
                %subjectdconn = cifti_conn_matrix(subject_dt_series,series,motion, FD_threshold, TR, minutes_limit, smoothing_kernal,L_surface,R_surface,bit8, remove_outliers, additional_mask);
                %temp_name = cifti_conn_matrix(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file, bit8, re')
                cii  = ciftiopen(subject_dt_series,wb_command);
                dtseries = cii.cdata;
                num_greys = size(dtseries,1);

                if num_greys == 91282
                    disp(['The number of greyordinates is ' num2str(num_greys) '. It is highly likely that there are subcortical voxels.  Removing...'])
                    disp('Opening cortex-only cifti for writing.')
                    cii  = ciftiopen(settings.path{11},wb_command); %path to surface_only dtseries.
                    dtseries=dtseries(1:59412,:);
                    cii.cdata = dtseries;
                    
                    ciftisave(cii, [output_directory filesep filename_short '_surf_only.dtseries.nii'], wb_command);
                    %ciftisave(cii, [path_to_orig_dtseries filesep filename_short '_surf_only.dtseries.nii'], wb_command)
                    
                    subject_dt_series = [output_directory filesep filename_short '_surf_only.dtseries.nii'];
                    %subject_dt_series = [path_to_orig_dtseries filesep filename_short '_surf_only.dtseries.nii'];
                    
                            motion_exten = strsplit(motion, '.');
                            motion_exten = char(motion_exten(end));
                            other_motion_mask = ~strcmp('mat', motion_exten);
                            
                            if other_motion_mask ==0
                                cmd =(['cp ' motion ' ' output_directory filesep motion_name_only 'surf_only.mat']);
                                disp([cmd '. This file is only copied so that an existing saved mask will not be overwritten.'])
                                system(cmd);
                                motion = [output_directory filesep motion_name_only 'surf_only.mat'];
                            else
                                cmd =(['cp ' motion ' ' output_directory filesep motion_name_only 'surf_only.txt']);
                                disp([cmd '. This file is only copied so that an existing saved mask will not be overwritten.'])
                                system(cmd);
                                motion = [output_directory filesep motion_name_only 'surf_only.txt'];
                            end
                            
                elseif num_greys == 59412
                    disp(['The number of greyordinates is ' num2str(num_greys) '. It is highly likely that subcortical have already been removed.  No subsectioning necessary'])
                    
                else
                    disp(['The number of greyordinates is ' num2str(num_greys) '. This is neither the expected 91282, nor 59412.  Unclear how to section your data. Exiting...'])
                    return
                end
                
                %                     cmd = [wb_command ' -cifti-separate ' subject_dt_series ' COLUMN ' ' -metric CORTEX_LEFT ' filename_short '.L.func.gii  -metric CORTEX_RIGHT ' filename_short '.R.func.gii'];
                %                     system(cmd)
                %                     cmd = [wb_command ' -cifti-create-dense-timeseries ' filename_short '.dtseries.nii -left-metric ' filename_short '.L.func.gii -right-metric ' filename_short 'R.func.gii  -timestep ' FD_threshold];
                %                     system(cmd)
                %                 large_subjectdconn = subjectdconn;
                %                 [subjectdconn] = surface_only_dconn(subjectdconn,'inferred');
                %                 disp('Removing connectivity matrix  since a smaller on eiwth surface only has been saved.')
                %cmd = (['rm -f ' num2str(large_subjectdconn)]);
                %system(cmd)
                already_surface_only =1;
            end
        end
        
        %% BUILD DCONN
        %subjectdconn = cifti_conn_matrix(subject_dt_series,series,motion, FD_threshold, TR, minutes_limit, smoothing_kernal,L_surface,R_surface,bit8, remove_outliers, additional_mask);
        subjectdconn = cifti_conn_matrix_for_wrapper_continous(wb_command, subject_dt_series, series, motion, FD_threshold, TR, minutes_limit,smoothing_kernal, left_surface_file, right_surface_file, bit8, remove_outliers, additional_mask, make_dconn_conc, [output_directory filesep], dtseries_conc, use_continous_minutes, memory_limit_value);
        %temp_name = cifti_conn_matrix(dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file, bit8, remove_outliers, additional_mask)
        %temp_name = cifti_conn_matrix   (dt_or_ptseries_conc_file,series,motion_file, FD_threshold, TR, minutes_limit, smoothing_kernal,left_surface_file, right_surface_file, bit8, remove_outliers, additional_mask)
        [path_to_dconn, dconn_name] = fileparts(subjectdconn);
        
        %% RUN COMMUNITY DECTECTION
        if run_infomap_too ==1
            
            %generate POTENTIAL filenames that might be removed later.
            file_split = strsplit(subjectdconn,'/');
            file_name = char(file_split(end));
            clu_split = strsplit(file_name,'.dc');
            cluname = char(clu_split(1));
            cluname = [cluname 'surface_only'];
            subjectdconn_surf_only =[cluname '.dconn.nii'];
            
            if exist([expected_final_dtseries_file],'file') == 2
                disp(['Final infomap dtseries file found: ' expected_final_dtseries_file])
                disp('Skipping infomap step.')
            else
                % some hardcodes for infomap
                distance_matrix=settings.path{19}; % path to distance matrix
                %some hardcodes for infomap:
                template = 'none';
                min_distance = 20;
                tie_density_vec = 'short';
                min_network_size = 400;
                min_region_size = 30;
                %community_detection = 'infomap';
                num_reps = 20;
                if strcmp(transform_data,'Convert_to_Zscores') ==1
                    donotZscore = 0;
                else
                    donotZscore = 1;
                end
                %surface_only = 0;
                %addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair'))
                addpath(genpath('/home/faird/shared/code/internal/utilities/community_detection/fair'))               
                [cleaned_info_ties_dtseries, raw_info_ties_dtseries] = Run_infomap_at_many_densities(subjectdconn,distance_matrix,template,min_distance,tie_density_vec,min_network_size,min_region_size,donotZscore,surface_only,already_surface_only,'infomap',num_reps);
                
               
                if surface_only ==1 % This name was likely altered during infomap due to the way it tries to remove the surface.
                    Comm_det_subject_folder = subjectdconn_surf_only;
                else
                    Comm_det_subject_folder = subjectdconn;
                end
                
                %            cmd=['ls -1d `pwd`' filesep 'Community_Detection_Min_Dist_20_TieDen_0.0*_MinNet_Size_400_MinReg_Size_30' filesep subjectdconn filesep 'wb_ready_files' filesep subID '_*_' num2str(smoothing_kernal) '.dtseries.nii_*_at_FD_' num2str(FD_threshold) '*perc_' num2str(min_distance) 'dist_bin_No_COMMUNITIES*dtseries.nii | sort '];
                %                [status, linux_results] = system(cmd,'-echo');
                %                 if status ==0
                %                     all_ties_singlecell = textscan(linux_results, '%s', 'delimiter', '\n' );
                %                     all_ties = all_ties_singlecell{1};
                %                 else
                %                     disp('Error in getting list of dtseries for each tie density.')
                %                     return
                %                 end
                
                for m=1:length(raw_info_ties_dtseries)
                    if exist(raw_info_ties_dtseries{m}) == 0
                        disp(['NOTE: Subject infomap series ' num2str(m) ' does not exist']);
                        return
                    else
                    end
                end
                
                disp('***MERGING ALL TIE DENISTIES INTO DTSERIES****')
                if use_all_ABCD_tasks == 1
                    MergeTimeSeries('TimeSeriesFiles',raw_info_ties_dtseries,'OutputFile',[cifti_output_folder filesep subID '_merged_tasks_infomap_densities.dtseries.nii'])
                    disp('***CALCULATING CONSENSUS FROM TIE DENISTIES AND CLEANING UP TINY PIECES****')
                    % Get concensus across densities and clean up little bits.
                    [merged_infomap_dtseries_vector] = clean_dscalars_by_size([cifti_output_folder filesep subID '_merged_tasks_infomap_densities.dtseries.nii'],[],[],[],[],30,[],1,1,0);
                else
                    MergeTimeSeries('TimeSeriesFiles',raw_info_ties_dtseries,'OutputFile',[cifti_output_folder filesep subID '_infomap_densities.dtseries.nii'])
                    disp('***CALCULATING CONSENSUS FROM TIE DENISTIES AND CLEANING UP TINY PIECES****')
                    % Get concensus across densities and clean up little bits.
                    [merged_infomap_dtseries_vector] = clean_dscalars_by_size([cifti_output_folder filesep subID '_infomap_densities.dtseries.nii'],[],[],[],[],30,[],1,1,0);
                    
                end
                
                %merged_infomap_dtseries_vector = [merged_infomap_dtseries_vector '.dtseries.nii'];
                disp('Final merged file has been sucesffully created. Removing  intermediate infomap files.')
                cmd = ['rm -fr `pwd`' filesep 'Community_Detection_Min_Dist_20_TieDen_0.0*_MinNet_Size_400_MinReg_Size_30' filesep Comm_det_subject_folder];
                disp(cmd);
                system(cmd)
            end
        else
            disp('Parameters are set to skip infomap')
        end
        
        %Step 2: get network assingments via template matching.
        [ eta_net_assign{i}, output_cifti_scalar_name] = template_matching_RH(subjectdconn, data_type, template_path,transform_data,output_cifti_name, cifti_output_folder ,wb_command,make_cifti_from_results, allow_overlap,overlap_method,surface_only,already_surface_only);
        
        
        if clean_up_intermed_files ==1 % RH added in case filespace becomes an limited.
            % Step 2.5: Remove dconn to save space.
            cmd = ['rm -f ' subjectdconn];
            disp(cmd);
            system(cmd);
            if use_all_ABCD_tasks == 1
                cmd = ['rm -f ' output_directory filesep filename_short '.dtseries.nii'];
                disp([cmd ' . Removing merged time series that created as part of the twins_mapping_wrapper.']);
                system(cmd);
                
                cmd = ['rm -f ' motion_path filesep motion_name_only '.mat'];
                
                disp([cmd ' . Removing motion mask that was created as part of the twins_mapping_wrapper.']);
                system(cmd);
            else
                disp('Keeping dconn. Be mindful of space.')
            end
            
            if surface_only ==1 && already_surface_only ==1 && clean_up_intermed_files == 1 % remove surface_only timeseries.
                % Step 2.5: Remove dconn to save space.
                if run_infomap_too ==1
                    cmd = ['rm -f ' [path_to_dconn filesep subjectdconn_surf_only]];
                    disp([cmd ' . This file will only be removed if it was created as part of the twins_mapping_wrapper.']);
                    system(cmd);
                end
                
                cmd = ['rm -f ' output_directory filesep filename_short '_surf_only.dtseries.nii'];
                disp([cmd ' . This file will only be removed if it was created as part of the twins_mapping_wrapper.']);
                system(cmd);
                
                cmd = ['rm -f ' output_directory filesep filename_short '_surf_only_SMOOTHED_' num2str(smoothing_kernal) '.dtseries.nii'];
                disp([cmd ' . This file will only be removed if it was created as part of the twins_mapping_wrapper.']);
                system(cmd);
              if other_motion_mask ==0  
                cmd = ['rm -f ' output_directory filesep motion_name_only 'surf_only.mat'];          
                disp([cmd ' . This file will only be removed if it was created as part of the twins_mapping_wrapper.']);
                system(cmd);
              else
              end
            else
                disp('Surface only files are not set to be removed. (Reduced time series and motion file).  You probably started with surface_only files, but check to make sure that excess files not being left behind.')
            end
            % Save motion mask
            if isnumeric(minutes_limit) ==1
                if surface_only ==1
                    cmd = ['mv -v ' motion_path filesep motion_name_only 'surf_only.mat_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt ' cifti_output_folder '_motion_masks/' ];
                    disp([cmd ' . Saving motion mask that was created as part of the twins_mapping_wrapper.']);
                    system(cmd);
                    
                else
                    cmd = ['mv -v ' motion_path filesep motion_name_only '.mat_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt ' cifti_output_folder '_motion_masks/' ];
                    disp([cmd ' . Saving motion mask that was created as part of the twins_mapping_wrapper.']);
                    system(cmd);
                end
            else % 'all frames was probably selected selected'
                if surface_only ==1
                    cmd = ['mv -v ' motion_path filesep motion_name_only 'surf_only.mat_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt ' cifti_output_folder '_motion_masks/' ];
                    disp([cmd ' . Saving motion mask that was created as part of the twins_mapping_wrapper.']);
                    system(cmd);
                    
                else
                    cmd = ['mv -v ' motion_path filesep motion_name_only '.mat_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt ' cifti_output_folder '_motion_masks/' ];
                    disp([cmd ' . Saving motion mask that was created as part of the twins_mapping_wrapper.']);
                    system(cmd);
                end
            end
            %sub-NDARINV3MTP07E9_ses-baselineYear1Arm1_merged_tasks_motionsurf_only.mat_0.2_cifti_censor_FD_vector_All_Good_Frames.txt
        end
    end %finish running template matching.
    
    %% Make pretty pictures of your results
    disp('Make pretty pictures of your results')
    if run_infomap_too == 1
        if use_all_ABCD_tasks == 1
            expected_final_dtseries_file = [cifti_output_folder filesep subID '_merged_tasks_infomap_densities_recolored.dtseries.nii'];
            if exist(expected_final_dtseries_file,'file') == 2
                merged_infomap_dtseries_vector=expected_final_dtseries_file;
            else
                expected_info_dscalar_name = [expected_final_dtseries_file(1:end-13),'.dscalar.nii'];
                if exist(expected_info_dscalar_name,'file') == 2
                    disp('Infomap dscalar file found.')
                else
                    disp(['Expected final infomap dtseries file not found. Expected to find: ' expected_final_dtseries_file])
                    disp('Since the TM results have been found and no infomap results have been found (the TM communites are calculated after running infomap), You may need to delete the template matching results and rerun.')
                    return
                end
                
            end
        else
            expected_final_dtseries_file = [cifti_output_folder filesep subID '_infomap_densities_recolored.dtseries.nii'];
            if exist(expected_final_dtseries_file,'file') == 2
                merged_infomap_dtseries_vector=expected_final_dtseries_file;
            else
                expected_info_dscalar_name = [expected_final_dtseries_file(1:end-13),'.dscalar.nii'];
                if exist(expected_info_dscalar_name,'file') == 2
                    disp('Infomap dscalar file found.')
                else
                    disp(['Expected final infomap dtseries file not found. Expected to find: ' expected_final_dtseries_file])
                    disp('Since the TM results have been found and no infomap results have been found (the TM communites are calculated after running infomap), You may need to delete the template matching results and rerun.')
                    return
                end
            end
        end
        
        disp('Converting final infomap dtseries to dscalar.')
        expected_info_dscalar_name = [expected_final_dtseries_file(1:end-13),'.dscalar.nii'];
        if exist(expected_info_dscalar_name,'file') == 2
            disp('Infomap dscalar file found. No conversion necessary')
        else
            if strcmp(merged_infomap_dtseries_vector(end-12:end),'.dtseries.nii') ==1 % take off .dtseries.nii extension for dscalar conversion
                merged_infomap_dtseries_vector = merged_infomap_dtseries_vector(1:end-13);
            else
                % The extension was probably removed when infomap was run, but if there
                % was error and you're rerunning the code, the extension might
                % still be on the file name.
            end
            cii_dscalar=ciftiopen([merged_infomap_dtseries_vector '.dtseries.nii'],wb_command);
            %newcii=cii_dscalar.cdata;
            
            %ciftisave(newcii,output_file,path_wb_c); % Making your cifti
            % modified from save to include reset-scalars flag.
            save(cii_dscalar,[merged_infomap_dtseries_vector '.gii'],'ExternalFileBinary')
            system([wb_command ' -cifti-convert -from-gifti-ext ' merged_infomap_dtseries_vector '.gii ' merged_infomap_dtseries_vector '.dscalar.nii -reset-scalars' ]);
            system([' rm -f ' merged_infomap_dtseries_vector '.gii ' merged_infomap_dtseries_vector '.dat ']);
            system([' rm -f ' merged_infomap_dtseries_vector '.dtseries.nii']);
            disp(['Done, file ' merged_infomap_dtseries_vector '  saved']);
            
        end
    end
    % Pictures defaults
    % The setting are hardcoded to make a power_colors png with 6 views. at 118 dots
    % per cm.
    
    %pics_code_path = '/home/faird/shared/code/internal/utilities/figure_maker/make_dscalar_pics_v9.3.sh';
    pics_code_path = settings.path{15}; % path to figure_maker bash script.
    pics_folder = [cifti_output_folder '/pics_template_matching'];
    disp(['mkdir -p ' pics_folder])
    system(['mkdir -p ' pics_folder]);
    if surface_only ==1
        make_subcortical_images = 'FALSE';
    else
        make_subcortical_images = 'TRUE';
    end
    
    %make template matching pic
    cmd = [pics_code_path ' ' recolored_name ' ' output_cifti_name '_recolored_TM ' pics_folder ' FALSE 1 18 power_surf FALSE 0 20 THRESHOLD_TEST_SHOW_OUTSIDE TRUE  ' make_subcortical_images ' png 8 118 FALSE ' wb_command ' ' settings.path{13} ' ' settings.path{14}];
    disp(cmd);
    system(cmd);
    
    %make infomap pic
    if run_infomap_too ==1
        pics_folder = [cifti_output_folder '/pics_infomap'];
        system(['mkdir -p ' pics_folder])
        cmd = [pics_code_path ' ' expected_info_dscalar_name ' ' output_cifti_name_info '_recolored_infomap ' pics_folder ' FALSE 1 18 power_surf FALSE 0 20 THRESHOLD_TEST_SHOW_OUTSIDE TRUE  ' make_subcortical_images ' png 8 118 FALSE ' wb_command ' ' settings.path{13} ' ' settings.path{14}];
        disp(cmd);
        system(cmd);
    end
    
    %%
    if calculate_mutual_info == 1
        disp('Done calculating network assingments for all pairs. Proceeding to calculate mutual information for all pairs.');
        
        for i = 1:2:num_sub %number of subjects
            %Step 3: calculate mutualinformation
            %idx = 1:2:26;
            %muI(i,1) = MutualInformation(eta_net_assign{i},eta_net_assign{i+1}); %Mutual information
            %[VIn(i,1), MIn(i,1)] = partition_distance(eta_net_assign{i},eta_net_assign{i+1}); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            muI(round(i/2),1) = MutualInformation(eta_net_assign{i},eta_net_assign{i+1}); %Mutual information
            [VIn(round(i/2),1), MIn(round(i/2),1)] = partition_distance(eta_net_assign{i},eta_net_assign{i+1}); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
        end
        
        %display summary
        disp(MIn)
        disp(muI)
        disp(VIn)
        
        %dlmwrite([cifti_output_folder '/' output_cifti_name '.txt'],MIn,muI,VIn,'roffset',1)
        save([cifti_output_folder '/' output_file_name '.mat'],'MIn','muI','VIn');
        disp('done calculating mutual information for all pairs')
    end
    
end
disp('Template_matching completed.')

end
