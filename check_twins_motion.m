function [all_subjects_minutes,all_mean_FD] = check_twins_motion(all_motion_conc_file, TR_or_TR_conc, FD_column,FD_conc_file,get_mean_FD,use_outlierdetection_mask_if_possible)

%This function read either a txt file (0 t d iscard , 1 to keep) or a Power_2014_FD_only.mat file.
% It will then tell you how many frames are available.
% TR te repitition time of your data.  For ADHD data set use 2.5
% THe FD column to use from the power_2014_FD-only file. if txt is used.The FD column won't be used.

%motion_all = importdata('/home/exacloud/lustre1/fnl_lab/projects/VA_Studies/Seed_map_3_studies/motion_concfile/MARC_MIND_PHN_motion.conc');
%motion_all = importdata('/home/groups/brainmri/infant/ECHO/PITT/unprocessed/niftis/sub-CT00434-2/ses-20190811/func/sub-CT00434-2_ses-20190811_task-rest_run-01_bold.json');
%all_motion_conc = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/split_halves/half1/Allframescensorhalf1.conc');

%Harcode warning
split_motion =1;
splits = 3;
save_split_motion_vectors = 1;
get_FDnumbers =1;

if iscell(all_motion_conc_file) ==1
    conc = 'false';
else
    conc = strsplit(all_motion_conc_file, '.');
    conc = char(conc(end));
end
if strcmp('conc',conc) == 1
    all_motion_conc = importdata(all_motion_conc_file);
else
    if iscell(all_motion_conc_file) ==0
        all_motion_conc = {all_motion_conc_file};
    else
        all_motion_conc = all_motion_conc_file;
    end
end
if exist('FD_conc_file', 'var') == 1 && ~isempty(FD_conc_file) == 1
    if strcmp('conc',conc) == 1
        FD_conc_file = importdata(FD_conc_file);
    else
        FD_conc_file = {FD_conc_file};
    end
else
end

if isnumeric(TR_or_TR_conc) == 1 % TR is one number (e.g. 2.5)
    if size(TR_or_TR_conc,1)>1
        TR_list = TR_or_TR_conc;
        many_TRs =1;
    else
        TR = TR_or_TR_conc;
        many_TRs =0;
    end
else
    if exist(TR_or_TR_conc,'file') == 1 %assume all the TRs are in a file.
        TR_list = importdata(TR_or_TR_conc);
        many_TRs =1;
    else
        TR = num2str(TR_or_TR_conc); % TR is a string
        many_TRs =0;
    end
end


%% Read motion file
%The number 31 here means FD =0.3, for FD=0.2 use 21
%TR = 2.2;
[~,~,c] = fileparts(all_motion_conc{1});
tic;
    good_frames = zeros(length(all_motion_conc),1);
    all_mean_FD = zeros(length(all_motion_conc),1);
    
for i = 1:length(all_motion_conc)
    
    if rem(i,100)==0
        disp([' Loading subject ' num2str(i)]);toc;
    end

    
    if strcmp(c,'.txt')
        FD_vector = importdata(all_motion_conc{i});
        good_frames(i) = sum(FD_vector);
        all_mean_FD(i,1) = 0;
    else
        
        
        load(all_motion_conc{i})
        if use_outlierdetection_mask_if_possible ==1
            try
                FD02_vector = motion_data{FD_column}.combined_removal;
            catch
                disp('Subject does not have the outlier dection mask in this .mat file. combined_removal mask not found. using power2014 FD mask instead.')
                FD02_vector =  motion_data{FD_column}.frame_removal; % use 21 for FD= 0.2, 31 for 0.3
            end
        end
        good_frames(i) = sum(abs(FD02_vector-1));
        if exist('get_mean_FD', 'var') == 1 && ~isempty(get_mean_FD) == 1
            all_mean_FD(i,1) = motion_data{1,FD_column}.remaining_frame_mean_FD;
            
        end
        disp(i);
    end
end

disp('Number of good frames;')
disp(good_frames);

disp('Amount of good data in minutes')

%% Calculate available minutes.
    all_subjects_minutes = zeros(length(all_motion_conc),1);
for i = 1:length(all_motion_conc)
    if strcmp(c,'.txt') % if TR is available from .mat file, use it to calculate the minutes of good data.
        if many_TRs ==0
            good_frames_in_minutes=good_frames(i,1)*TR/60;
            poss_frames = length(importdata(all_motion_conc{i}));
            possible_minutes = poss_frames*TR/60;
        else
            good_frames_in_minutes=(good_frames(i,1)*TR_list(i))/60;
            poss_frames = length(importdata(all_motion_conc{i}));
            possible_minutes = (poss_frames*TR_list(i))/60;
        end
        
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
        all_subjects_minutes(i,1) = good_frames_in_minutes;
    else
        load(all_motion_conc{i})
        good_frames_in_minutes=good_frames(i,1)*motion_data{1,FD_column}.epi_TR/60;
        poss_frames =  length(motion_data{FD_column}.frame_removal);
        possible_minutes = poss_frames*motion_data{1,FD_column}.epi_TR/60;
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
        all_subjects_minutes(i,1) = good_frames_in_minutes;
        
        
    end
end



%%Calculate pre  and post scrubbing  average FD
if exist('FD_conc_file', 'var') == 1 && ~isempty(FD_conc_file) == 1
    
    if exist('FD_conc_file', 'var') == 1 && FD_conc_file == 1
        for i = 1:length(all_motion_conc)
            load(FD_conc_file{i})
            FD = motion_numbers.FD;
            mask = importdata(all_motion_conc{i});
            FD_afterscrub = FD(logical(mask));
            avgFD_prescrub(i,:) =mean(FD);
            avgFD_postscrub(i,:) = mean(FD_afterscrub);
        end
    else
        
    end
    
    %%Run data
    if exist('find_best_subjects', 'var') == 1 && find_best_subjects == 1
        [frames, best_subs] = sort(good_frames');
        
        upper_threshold = 346; % %provide a list of subjects that have frames above this threshold. e.g. 346 frames = ~14.4 minutes of data.
        thresidx = find(frames ==upper_threshold);
        thresidx = thresidx(1); % if many subjets are exactly at the threshold, take the first.
        for i = thresidx:length(frames)
            lots_o_data_subs{i-thresidx+1} = all_motion_conc{best_subs(i)};
        end
        
        
        j = 1;
        if exist('split_motion', 'var') == 1 && split_motion == 1
            for i = 1:length(lots_o_data_subs)
                load(lots_o_data_subs{i})
                FD_vector =  motion_data{31}.frame_removal;
                for k = 1:splits %3 splits
                    restsize = round(length(FD_vector)/splits);
                    Rest(:,k)=FD_vector((restsize*k-restsize+1):(restsize*k));
                    if exist('save_split_motion_vectors,var') == 1 && save_split_motion_vectors ==1
                    else
                        %save a list of files being written
                        allrestlist{j} = [fileparts(fileparts(fileparts(lots_o_data_subs{i}))) '/MNINonLinear/Results/REST' num2str(k) '/FNL_preproc/REST_split_FD_mask.txt'];j=j+1;
                        %save files
                        
                        %dlmwrite([fileparts(fileparts(fileparts(lots_o_data_subs{i}))) '/MNINonLinear/Results/REST' num2str(k) '/FNL_preproc/REST_split_FD_mask.txt'],Rest(:,k))
                    end
                end
            end
            %disp('Saving list of motion vectors...')
            %dlmwrite('RESTfiles.txt',allrestlist')
        else
        end
    else
    end
else
end
disp('Done')
