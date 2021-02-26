function [tasks, tasks_motions] = make_scan_conc(input_directory, timeseriesname)
%This code works by taking the filename and time series of a particular
%scan (e.g. the task-rest) and finding the other tasks.
%03/04/2020 -Robert Hermosillo

%Currently the 4 tasks to look for the REST , MID, SST, nback tasks.


z =1; %Counter for found tim series
x=1; %Counter for found motion

%% Check for rests

disp(['Checking for all possible task/rest scans in the following directory: ' input_directory])
[status, task_rest_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_bold_timeseries.dtseries.nii']);
if status ==0
    disp('rest time series found')
    file_singlecell = textscan(task_rest_dtseries, '%s', 'delimiter', '\n' );
    task_rest_dtseries = file_singlecell{1};
    tasks{z,1} = task_rest_dtseries;    z= z+1;
else
    disp('No resting state scan found')
end
%check for corresponding motion file
[status, task_rest_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-rest_bold_mask.mat']);
if status ==0
    disp('rest motion found')
    file_singlecell = textscan(task_rest_motion, '%s', 'delimiter', '\n' );
    task_rest_motion = file_singlecell{1};
    tasks_motions{x,1} = task_rest_motion;     x= x+1;
else
    disp('No resting motion .mat found')
end

%% Check for MID
[status, task_MID_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-MID_bold_timeseries.dtseries.nii']);
if status ==0
    disp('MID time series found')
    file_singlecell = textscan(task_MID_dtseries, '%s', 'delimiter', '\n' );
    task_MID_dtseries = file_singlecell{1};
    tasks{z,1} = task_MID_dtseries;    z= z+1;
else
    disp('No MID time series scan found')
end
%check for corresponding motion file
[status, task_MID_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-MID_bold_mask.mat']);
if status ==0
    disp('MID motion found')
    file_singlecell = textscan(task_MID_motion, '%s', 'delimiter', '\n' );
    task_MID_motion = file_singlecell{1};
    tasks_motions{x,1} = task_MID_motion;     x= x+1;
else
    disp('No MID motion .mat found')
end


%% Check fot SST
[status, task_SST_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-SST_bold_timeseries.dtseries.nii']);
if status ==0
    disp('SST time series found')
    file_singlecell = textscan(task_SST_dtseries, '%s', 'delimiter', '\n' );
    task_SST_dtseries = file_singlecell{1};
    tasks{z,1} = task_SST_dtseries;    z= z+1;
else
    disp('No SST time series scan found')
end
%check for corresponding motion file
[status, task_SST_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-SST_bold_mask.mat']);
if status ==0
    disp('SST motion found')
    file_singlecell = textscan(task_SST_motion, '%s', 'delimiter', '\n' );
    task_SST_motion = file_singlecell{1};
    tasks_motions{x,1} = task_SST_motion;     x= x+1;
else
    disp('No SST motion .mat found')
end


%% check for nBack
[status, task_nback_dtseries] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_bold_timeseries.dtseries.nii']);
if status ==0
    disp('nback time series found')
    file_singlecell = textscan(task_nback_dtseries, '%s', 'delimiter', '\n' );
    task_nback_dtseries = file_singlecell{1};
    tasks{z,1} = task_nback_dtseries;    z= z+1;
else
    disp('No nback time series scan found')
end
%check for corresponding motion file
[status, task_nback_motion] = system(['ls ' input_directory filesep '*_ses-baselineYear1Arm1_task-nback_bold_mask.mat']);
if status ==0
    disp('nback motion found')
    file_singlecell = textscan(task_nback_motion, '%s', 'delimiter', '\n' );
    task_nback_motion = file_singlecell{1};
    tasks_motions{x,1} = task_nback_motion;     x= x+1;
else
    disp('No nback motion .mat found')
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