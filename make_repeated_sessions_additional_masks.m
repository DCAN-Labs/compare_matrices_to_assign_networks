function [all_mask_vec] = make_repeated_sessions_additional_masks(input_directories, wb_command, half_to_use,check_mask_against_merged_dtseries)

%This code checks for run-level dtseries. For run numerber that are larger than 10, and makes a motion
%masks of runs that are smaller or larger than 10.  It will then write out
%a mask (of 1s and 0s) that correspond with runs greater than or less than
%0.
%03/04/2020 -Robert Hermosillo

%inputs are:
%1) input_directory (char): the path to input directory that contains the dtseries(s) files that you want to use.
%2) combine_all_dtseries_in_dir (numeric 1 or 0):if 1 the code will use all dtseries in the current folder ending in "*desc-filtered_timeseries.dtseries.nii".
%3) split_runs (numeric 1 or 0): if set to 1, the code will automatically split runs with runs 1-10 and 11-20 into
%4) halfmasktoexport (1 or 2) :the code generate will generate a mask of 1s and zeros
%based on run numbers and export that mask to use later.

%If combine_all_dtseries_in_dir is set to 1, the code will look for all dtseries.  If set to 0, it will look for the following ABCD tasks:
%Currently the 4 tasks to look for the REST , MID, SST, nback tasks.

task_root_list = {'gngreg';'gngrew';'rest'};

if half_to_use ==1
elseif half_to_use ==2
else
    error('half_to_use must be either 1 or 2.')
end

conc = strsplit(input_directories, '.');
conc = char(conc(end));
if strcmp('conc', conc)
    input_directories = importdata(input_directories);
else
    input_directories = {input_directories};
end


for d =1:size(input_directories)
    input_directory= char(input_directories(d));
    
    if check_mask_against_merged_dtseries ==1
        [status, merged_all_task_dtseries_half1] = system(['ls ' input_directory filesep '*merged_tasks.dtseries.nii']);
        if status ~= 0
            error(['Check dtseries file: ' merged_all_task_dtseries_half1])
        end
        cmd = [wb_command ' -file-information -only-number-of-maps '  merged_all_task_dtseries_half1];
        disp(cmd)
        [~,length_merged_tseries] = system(cmd);
        length_merged_tseries = str2double(length_merged_tseries);
    end
    
    for i = 1:size(task_root_list,1)
        task_root = char(task_root_list(i));
        
        [status, all_task_dtseries_half1] = system(['ls ' input_directory filesep '*' task_root '*run-0*_bold_timeseries.dtseries.nii']);
        
        if status ==0
            disp([' Time series found for task: ' char(task_root_list(i))]);
            file_singlecell1 = textscan(all_task_dtseries_half1, '%s', 'delimiter', '\n' );
            num_dtseries1 = length(file_singlecell1{1,1});
            disp(['number of dtseries found (half1): ' num2str(num_dtseries1)]);
            %task_rest_dtseries = file_singlecell{1};
        else
            disp(['No resting state scan found for task: ' char(task_root_list(i))])
            file_singlecell1 = {};
            num_dtseries1 = 0;
        end
        
        [status, all_task_dtseries_half2] = system(['ls ' input_directory filesep '*' task_root '*run-1*_bold_timeseries.dtseries.nii']);
        if status ==0
            disp([' Time series found for task:' char(task_root_list(i))]);
            file_singlecell10 = textscan(all_task_dtseries_half2, '%s', 'delimiter', '\n' );
            num_dtseries2 = length(file_singlecell10{1,1});
            disp(['number of dtseries found: ' num2str(num_dtseries2)]);
            %task_rest_dtseries = file_singlecell{1};
        else
            disp(['No resting state scan found for task (half2): ' char(task_root_list(i))])
            file_singlecell10 = {};
            num_dtseries2 = 0;
        end
        
        
        if num_dtseries1 ==0
            this_task_runs01 =[];
        else
            for j=1:num_dtseries1
                cmd = [wb_command ' -file-information -only-number-of-maps '  file_singlecell1{1,1}{j,1}];
                disp(cmd)
                [~,length_tseries] = system(cmd);
                disp(char(length_tseries))
                length_tseries = str2double(length_tseries);
                
                if length_tseries <=6
                    disp(['WARNING: ' char(file_singlecell1{1,1}{j,1}) ' tseries is short and will not be used.']);
                    length_tseries =0;
                end
                if j ==1
                    this_task_runs01 = ones(length_tseries,1);
                else
                    new_task_runs1 = ones(length_tseries,1);
                    this_task_runs01 = [this_task_runs01; new_task_runs1];
                end
            end
        end
        if num_dtseries2 ==0
            this_task_runs10 =[];
        else
            for k=1:num_dtseries2
                cmd = [wb_command ' -file-information -only-number-of-maps '  file_singlecell10{1,1}{k,1}];
                disp(cmd)
                [~,length_tseries] = system(cmd);
                disp(num2str(length_tseries));
                length_tseries = str2double(length_tseries);
                
                if length_tseries <=6
                    disp(['WARNING: ' char(file_singlecell1{1,1}{k,1}) ' tseries is short and will not be used.']);
                    length_tseries =0;
                end
                if k ==1
                    this_task_runs10 = zeros(length_tseries,1);
                else
                    new_task_runs2 = zeros(length_tseries,1);
                    this_task_runs10 = [this_task_runs10; new_task_runs2];
                end
            end
        end
        
        %% Build mask for each task
        if i ==1
            task01_mask = [this_task_runs01; this_task_runs10];
        end
        if i==2
            task02_mask = [this_task_runs01; this_task_runs10];
        end
        if i==3
            task03_mask = [this_task_runs01; this_task_runs10];
        end
        
        allvectors01{i,1} =  this_task_runs01;
        allvectors10{i,1} = this_task_runs10;
        
        clear file_singlecell1 file_singlecell10 this_task_runs01 this_task_runs10% clear this cell so that previously founc files are not incorrectly assumed to be found.
    end
    
    
    all_mask_vec = [task01_mask; task02_mask; task03_mask];
    if half_to_use ==1
    else
        all_mask_vec  = ~all_mask_vec ;
    end
    
    if  check_mask_against_merged_dtseries ==1
        if length(all_mask_vec) ~= length_merged_tseries
            disp(['dtseries length: ' num2str(length_merged_tseries)]);
            disp(['mask length: ' num2str(length(all_mask_vec))]);
            error(['the merged dtseries is not the same length of the runs. for sub: ' num2str(d) ]);
        else
            disp(['dtseries length: ' num2str(length_merged_tseries)]);
            disp(['mask length: ' num2str(length(all_mask_vec))]);
            disp('Writing mask')
            fileID = fopen([input_directory filesep 'merged_additional_mask_half' num2str(half_to_use) '.txt'],'w');
            fprintf(fileID,'%d\n',all_mask_vec);
            fclose(fileID);
        end
    else
        disp('Mask has not been checked.')
        disp(['mask length: ' num2str(length(all_mask_vec))]);
        disp('Writing mask')
        fileID = fopen([input_directory filesep 'merged_additional_mask_half' num2str(half_to_use) '.txt'],'w');
        fprintf(fileID,'%d\n',all_mask_vec);
        fclose(fileID);
    end
    disp(['mask is here: ' input_directory filesep 'merged_additional_mask_half' num2str(half_to_use) '.txt'])
end %input directories
end %function