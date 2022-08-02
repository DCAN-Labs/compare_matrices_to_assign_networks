function [tasks, tasks_motions, halfmasktoexport] = make_scan_conc(input_directory, combine_non_ABCD_dtseries_in_dir,split_runs,halfmasktoexport,varargin)

%This code works by taking the filename and time series of a particular
%scan (e.g. the task-rest) and finding the other tasks.
%03/04/2020 -Robert Hermosillo

%02/08/2022 - Robert Hermosillo

%Modified to look for all dtseries and motion within the supplied folder.
%inputs are:
%1) input_directory (char): the path to input directory that contains the dtseries(s) files that you want to use.
%2) combine_all_dtseries_in_dir (numeric 1 or 0):if 1 the code will use all dtseries in the current folder ending in "*desc-filtered_timeseries.dtseries.nii".
%3) split_runs (numeric 1 or 0): if set to 1, the code will automatically split runs with runs 1-10 and 11-20 into
%4) halfmasktoexport (1 or 2) :the code generate will generate a mask of 1s and zeros
%based on run numbers and export that mask to use later.
%varagin
if ~isempty(varargin)
    use_custom_pattern = varargin{1};
    ls_command_pattern = varargin{2};
else
    use_custom_pattern = 0;
    ls_command_pattern = '';
end
%If combine_all_dtseries_in_dir is set to 1, the code will look for all dtseries.  If set to 0, it will look for the following ABCD tasks:
%Currently the 4 tasks to look for the REST , MID, SST, nback tasks.

z =1; %Counter for found tim series
x=1; %Counter for found motion

%% Check for rests
if combine_non_ABCD_dtseries_in_dir ==1
    disp(['Checking for all possible task/rest scans in the following directory: ' input_directory])
    if use_custom_pattern ==1
        disp('Using custom string pattern to find dtseries.')
        [status, task_rest_dtseries] = system(['ls ' input_directory filesep ls_command_pattern]);
    else
        [status, concat_task_rest_dtseries] = system(['ls ' input_directory filesep '*desc-filtered_timeseries.dtseries.nii']);
        if status ==0
            disp('Concatenated Rest time series found')
            concat_file_singlecell = textscan(concat_task_rest_dtseries, '%s', 'delimiter', '\n' );
            concat_num_dtseries = length(concat_file_singlecell{1,1});
            disp(['number of concat dtseries found: ' num2str(concat_num_dtseries)]);
        else
            disp('No concatenated resting state scan found')
            
        end
    end
    
    if split_runs ==1
        
        for i=1: concat_num_dtseries
            task_split_name_cells = strsplit(concat_file_singlecell{1,1}{i,1},'_task');
            task_split_name_cells2 = strsplit(task_split_name_cells{1,1},'_');
            task_names{i,1} = task_split_name_cells2;
            [status, task_rest_dtseries] = system(['ls ' input_directory filesep task_name '*run*_timeseries.dtseries.nii']);
            
            if status ==0
                disp(' run timeseries found')
                file_singlecell = textscan(task_rest_dtseries, '%s', 'delimiter', '\n' );
                num_dtseries = length(file_singlecell{1,1});
                disp(['number of dtseries found: ' num2str(num_dtseries)]);
                halfmasktoexport=zeros(num_dtseries,1);
                %task_rest_dtseries = file_singlecell{1};
                for i=1:num_dtseries
                    this_dtserires_string = file_singlecell{1,1}{i,1};
                    my2filestrings=strsplit(this_dtserires_string,'_run-');
                    my2newfilestrings=strsplit(my2filestrings{2},'_bold_');
                    run_number=str2double(my2newfilestrings{1});
                    if half_to_export ==1
                        if run_number < 10 %assume that runs that are less than 10 are from the first half.
                            halfmasktoexport(i) =1;
                        else
                            halfmasktoexport(i) =0;
                        end
                    else %half is equal to 2
                        if run_number > 10 %assume that runs that are less than 10 are from the first half.
                            halfmasktoexport(i) =1;
                        else
                            halfmasktoexport(i) =0;
                        end
                    end
                    
                end
            else
                disp('No resting state scan found')
            end
        end
    end
    
    for i=1:concat_num_dtseries
        
        dtseries_name = concat_file_singlecell{1,1}{i,1};
        tasks{z,1} = dtseries_name;    z= z+1;
        [filedir,rootname_1ext] = fileparts(dtseries_name);
        [~,rootname_noext] = fileparts(rootname_1ext);
        
        disp('Assuming that dseries ends with "_timeseries".')
        short_root = char(rootname_noext(1:end-11));
        disp(['dtseries root name is ' char(short_root)])
        
        corresponding_motion_file = [short_root '_motion_mask.mat'];
        disp(['looking looking for .mat file that has the following name: ' filedir filesep corresponding_motion_file ])
        
        [status, task_motion] = system(['ls ' filedir filesep corresponding_motion_file]);
        if status ==0
            disp('rest motion found')
            motion_singlecell = textscan(task_motion, '%s', 'delimiter', '\n' );
            task_motion = motion_singlecell{1};
            tasks_motions{x,1} = task_motion;     x= x+1;
        else
            disp('No resting motion .mat found')
        end
    end
    
else
    disp(['Checking for ABCD possible task/rest scans in the following directory: ' input_directory])
    
    %[status, task_rest_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_bold_timeseries.dtseries.nii']);
     [status, task_rest_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_bold_desc-filtered_timeseries.dtseries.nii']);
   
    if status ==0
        disp('rest time series found')
        file_singlecell = textscan(task_rest_dtseries, '%s', 'delimiter', '\n' );
        task_rest_dtseries = file_singlecell{1};
        tasks{z,1} = task_rest_dtseries;    z= z+1;
    else
        disp('No resting state scan found')
    end
    %check for corresponding motion file
    %[status, task_rest_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_bold_mask.mat']);
     [status, task_rest_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_desc-filtered_motion_mask.mat']);
   
    if status ==0
        disp('rest motion found')
        file_singlecell = textscan(task_rest_motion, '%s', 'delimiter', '\n' );
        task_rest_motion = file_singlecell{1};
        tasks_motions{x,1} = task_rest_motion;     x= x+1;
    else
        disp('No resting motion .mat found')
    end
    
    %% Check for MID
    [status, task_MID_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-MID_bold_desc-filtered_timeseries.dtseries.nii']);
    if status ==0
        disp('MID time series found')
        file_singlecell = textscan(task_MID_dtseries, '%s', 'delimiter', '\n' );
        task_MID_dtseries = file_singlecell{1};
        tasks{z,1} = task_MID_dtseries;    z= z+1;
    else
        disp('No MID time series scan found')
    end
    %check for corresponding motion file
    %[status, task_MID_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-MID_bold_mask.mat']);
     [status, task_MID_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-MID_desc-filtered_motion_mask.mat']);
 
    if status ==0
        disp('MID motion found')
        file_singlecell = textscan(task_MID_motion, '%s', 'delimiter', '\n' );
        task_MID_motion = file_singlecell{1};
        tasks_motions{x,1} = task_MID_motion;     x= x+1;
    else
        disp('No MID motion .mat found')
    end
    
    
    %% Check fot SST
    [status, task_SST_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-SST_bold_desc-filtered_timeseries.dtseries.nii']);
    if status ==0
        disp('SST time series found')
        file_singlecell = textscan(task_SST_dtseries, '%s', 'delimiter', '\n' );
        task_SST_dtseries = file_singlecell{1};
        tasks{z,1} = task_SST_dtseries;    z= z+1;
    else
        disp('No SST time series scan found')
    end
    %check for corresponding motion file
    %[status, task_SST_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-SST_bold_mask.mat']);
     [status, task_SST_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-SST_desc-filtered_motion_mask.mat']);

    if status ==0
        disp('SST motion found')
        file_singlecell = textscan(task_SST_motion, '%s', 'delimiter', '\n' );
        task_SST_motion = file_singlecell{1};
        tasks_motions{x,1} = task_SST_motion;     x= x+1;
    else
        disp('No SST motion .mat found')
    end
    
    
    %% check for nBack
    %[status, task_nback_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_bold_desc-filtered_timeseries.dtseries.nii']);
    [status, task_nback_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_bold_desc-filtered_timeseries.dtseries.nii']);
 
    if status ==0
        disp('nback time series found')
        file_singlecell = textscan(task_nback_dtseries, '%s', 'delimiter', '\n' );
        task_nback_dtseries = file_singlecell{1};
        tasks{z,1} = task_nback_dtseries;    z= z+1;
    else
        disp('No nback time series scan found')
    end
    %check for corresponding motion file
    %[status, task_nback_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_bold_mask.mat']);
    [status, task_nback_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_desc-filtered_motion_mask.mat']);
  
    if status ==0
        disp('nback motion found')
        file_singlecell = textscan(task_nback_motion, '%s', 'delimiter', '\n' );
        task_nback_motion = file_singlecell{1};
        tasks_motions{x,1} = task_nback_motion;     x= x+1;
    else
        disp('No nback motion .mat found')
    end
end

disp(['Numer of tasks found: ' num2str(size(tasks,1))])
disp(['Numer of tasks found: ' num2str(size(tasks_motions,1))])
% check ls results
if size(tasks,1)  == size(tasks_motions,1)
    disp('Equal numbers of tasks and corresponding motion files found. This is a good thing.')
else
    disp('Different numbers of tasks and corresponding motion files found. Something is missing. Exiting...')
    return
end


end