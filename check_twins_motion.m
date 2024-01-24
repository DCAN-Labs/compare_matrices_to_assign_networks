function [all_subjects_minutes,all_mean_FD,allrestlist] = check_twins_motion(all_motion_conc_file, TR_or_TR_conc, output_dir, output_name, FD_column,FD_conc_file,get_mean_FD,use_outlierdetection_mask_if_possible,split_motion,splits)

%This function read either a txt file (0 t d iscard , 1 to keep) or a Power_2014_FD_only.mat file.
% It will then tell you how many frames (and minutes) are available per subject.
% If you additionally provide conc file of the FD values, it will calculate
% the pre-censor and post-censor average movement.

%inputs are: 
% all_motion_conc_file = A path to the power_2014_FD_only.mat or a .txt
% file of 1s and 0s(where 1 indicates keep and 0 indicates discard). (e.g. '/path/to/sub01_ses01_power2014_FDonly.mat' or '/path/to/all_subject_power2014_FDonly.conc'

% TR_or_TR_or_TR_conc = the BOLD repetition time (in seconds) (e.g. 0.8 or
% '/path/to/all_mysubjectTRs').

%output_dir = the path to the outputdirectory where you woudld like to save the data.   

%output_name = The name of the output file that you would like to use for
%the .mat file.

%FD_column = The FD column to use from the power_2014_FD-only file. if txt is used. The FD column variable won't be used.

%FD_conc_file = path to the FD conc file.  This will allow the code to
%calculate the amount of motion before and after motion censoring.  You can
%provide an emptry string if you don't have this file.

%get_mean_FD = set to 1, to get the mean FD, set to 0 if you don't want to
%calculate the mean FD.
%use_outlierdetection_mask_if_possible = in the new versions of dcan bold
%proc, outlier detection is calcuclated and saved as a vector.  You can use
%this mask.
%split_motion,splits)

%Harcode warning
%split_motion =1;
%splits = 2;
save_split_motion_vectors = 1;
get_FDnumbers =1;
find_subjects_with_most_data=1;
%minimum_frames_threshold = 750; % if find_subjects_with_most_data is ==1, then subset subjects with at least this many frames.
minimum_frames_threshold = 1500; % if find_subjects_with_most_data is ==1, then subset subjects with at least this many frames (@TR=0.8, 1500=20 minutes).

use_all_frames_in_splits = 1; % set to 1 if you want to use all available frames.  (eg. if have 6000 good frames, and you want to make 6 splits of 1000).  Setting this value to zero will make X splits where each is randomly sampled to fit the minimum threshold criteria (e.g you have 6000 frames and you want 2 splits of 80 frames.  In this instance you have way more available frames, so the code will cut data in half then randomly sample 80 frames.

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
        FD_vector_1isbad = importdata(all_motion_conc{i});
        good_frames(i) = sum(FD_vector_1isbad);
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
            possible_frames = length(importdata(all_motion_conc{i}));
            possible_minutes = possible_frames*TR/60;
        else
            good_frames_in_minutes=(good_frames(i,1)*TR_list(i))/60;
            possible_frames = length(importdata(all_motion_conc{i}));
            possible_minutes = (possible_frames*TR_list(i))/60;
        end
        
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
        all_subjects_minutes(i,1) = good_frames_in_minutes;
    else
        load(all_motion_conc{i})
        good_frames_in_minutes=good_frames(i,1)*motion_data{1,FD_column}.epi_TR/60;
        possible_frames =  length(motion_data{FD_column}.frame_removal);
        possible_minutes = possible_frames*motion_data{1,FD_column}.epi_TR/60;
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
        all_subjects_minutes(i,1) = good_frames_in_minutes;
        
        
    end
    %all_good_frames(i,1) = good_frames;
    all_possible_minutes(i,1) = possible_minutes;
    all_poss_frames(i,1) = possible_frames;
    
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
end

%Find subjects with the most data
if exist('find_subjects_with_most_data', 'var') == 1 && find_subjects_with_most_data == 1
    [frames, best_subs] = sort(good_frames);
    
    %upper_threshold = 346; % %provide a list of subjects that have frames above this threshold. e.g. 346 frames = ~14.4 minutes of data. 1500 at a TR of 0.8 = 1200sec. 1200secs / 60sec/min is 20minutes.
    %1500 at a TR of 0.8 = 1200sec. 1200secs / 60sec/min is 20minutes.
    %minimum_frames_threshold = 1208;
    thresidx = find(frames >=minimum_frames_threshold);
    if isempty(thresidx)
        thresidx = find(frames ==minimum_frames_threshold+1);
    end
    thresidx = thresidx(1); % if many subjects are exactly at the threshold, take the first.
    
    for i = thresidx:length(frames)
        lots_o_data_subs{i-thresidx+1,1} = all_motion_conc{best_subs(i)};
    end
    
    
    j = 1;
    %To generate a random mask with exactly 10 minutes in each split half, we first half to ...
    %1)find the mid point of good frames.
    %2)randomly sample the frames to make exactly x minutes (for each half)
    %3)build write those 1s to en empty mask.
    
    if exist('split_motion', 'var') == 1 && split_motion == 1
                 f1=figure;
         f2=figure;
         f3=figure;
        for i = 1:length(lots_o_data_subs)
            load(lots_o_data_subs{i})
            %FD_vector_1isbad =  motion_data{31}.frame_removal;
            
            if use_outlierdetection_mask_if_possible ==1
                try
                    FD_vector_1isbad = motion_data{FD_column}.combined_removal;
                    FD_vector_1isgood = double(~FD_vector_1isbad);
                catch
                    disp('Subject does not have the outlier dection mask in this .mat file. combined_removal mask not found. using power2014 FD mask instead.')
                    FD_vector_1isbad =  motion_data{FD_column}.frame_removal; % use 21 for FD= 0.2, 31 for 0.3
                    FD_vector_1isgood = double(~FD_vector_1isbad);
                     FD_vector_1isgood_idx = find(FD_vector_1isgood==1);
                   
                end
            end
            
            %restsize = round(length(FD_vector_1isgood)/splits);
            total_good_frames = sum(abs(FD_vector_1isbad-1));
            good_frames_per_split = floor(total_good_frames/splits); % round down to the nearst whole number.
            goodframes_sofar = cumsum(FD_vector_1isgood);
            Frame_cut = find(goodframes_sofar == good_frames_per_split);
            startframe_temp=1;
            
            for k = 1:splits % splits
                mask_vec = zeros(length(FD_vector_1isbad),1);
                if k ==1
                        startframe_temp=1;                  
                else
                    %startframe_temp = startframe_temp+(k-1)*Frame_cut;
                     startframe_temp = ((k-1)*good_frames_per_split)+1;
                   
                end
                %startframe_temp = startframe_temp+(k-1)*Frame_cut;
                if k ==1
                    %split_vec = FD_vector_1isgood(startframe_temp:Frame_cut*k);
                    %split_vec = FD_vector_1isgood(startframe_temp:frames_per_split*k);
                    split_vec = FD_vector_1isgood_idx(startframe_temp:good_frames_per_split*k);
                    
                else
                    %split_vec = FD_vector_1isgood(startframe_temp:end);
                    %split_vec = FD_vector_1isgood(startframe_temp:good_frames_per_split*k);
                    split_vec = FD_vector_1isgood_idx(startframe_temp:good_frames_per_split*k);
                    
                    
                end
                
                %good_frames_idx = find(FD_vector_1isbad==0);
                if k ==1
                    %mask_vec(startframe_temp:Frame_cut*k)=split_vec;
                     mask_vec(startframe_temp:good_frames_per_split*k)=split_vec;
                   
                else
                    mask_vec(startframe_temp:good_frames_per_split*k)=split_vec;
                end
                split_good_frames_idx = find(mask_vec ==1);
                
                
                %restsize = round(length(FD_vector)/splits);
                
                if use_all_frames_in_splits ==0
                    % Number of good frames to randomly pull
                    %rand_good_frames = sort(randperm(length(good_frames_idx), round(minutes_limit*60/TR)));
                    %rand_good_frames = sort(randperm(length(split_good_frames_idx), frames_per_split));
                    rand_good_frames = sort(randperm(length(split_good_frames_idx), floor(minimum_frames_threshold/splits)));
                    %FDvec_cut = zeros(length(FDvec), 1);
                    ones_idx = split_good_frames_idx(rand_good_frames);
                else
                    ones_idx =   split_good_frames_idx;
                end
                
                % The new vector that should match good frames with the
                % minutes limit
                mask_vec_min_lim = zeros(length(FD_vector_1isbad),1);
                mask_vec_min_lim(ones_idx) = 1;
                
                %mask_vec_min_lim(start_frame:start_frame+frames_per_split) = FDvec_cut;
                
                %Rest(:,k)=FD_vector_1isbad((restsize*k-restsize+1):(restsize*k));
                
                if exist('save_split_motion_vectors','var') == 1 && save_split_motion_vectors == 1
                    
                    %save a list of files being written
                    %allrestlist{j} = [fileparts(fileparts(fileparts(lots_o_data_subs{i}))) '/MNINonLinear/Results/REST' num2str(k) '/REST_split' num2str(k) '_FD_mask.txt'];
                    %j=j+1;
                    
                    figure(f1);
                    bar(FD_vector_1isgood)
                    figure(f2);
                    bar(mask_vec);
                    figure(f3);
                    bar(mask_vec_min_lim);
                    
                    
                    %save files
                    [~, file_root] =fileparts(lots_o_data_subs{i});
                    dlmwrite([output_dir filesep file_root '_split' num2str(k) '_FD' num2str(motion_data{1,FD_column}.FD_threshold) '_mask.txt'],mask_vec_min_lim)
                    allrestlist{j,1} = [output_dir filesep file_root '_split' num2str(k) '_FD' num2str(motion_data{1,FD_column}.FD_threshold) '_mask.txt'];
                    all_best_masks{j,1} = mask_vec_min_lim;
                    j=j+1;
                    clear split_good_frames_idx mask_vec mask_vec_min_lim ones_idx
                else
                end
                %startframe_temp = Frame_cut +1;
                %k=k+1;
            end
        end
        %disp('Saving list of motion vectors...')
        %dlmwrite('RESTfiles.txt',allrestlist')
    else
    end
    allrestlist = {};
    if exist('split_motion', 'var') == 1 && split_motion == 1
    save([output_dir filesep  output_name 'best_subs_split' num2str(splits) '_FD' num2str(motion_data{1,FD_column}.FD_threshold) '_masks.mat'],'allrestlist','all_best_masks','all_subjects_minutes','all_mean_FD','allrestlist','all_motion_conc_file')
    else
     save([output_dir filesep  output_name '_masks_summary.mat'],'minimum_frames_threshold','lots_o_data_subs','allrestlist', 'all_poss_frames', 'all_possible_minutes', 'all_subjects_minutes','all_mean_FD','allrestlist','all_motion_conc_file')       
    end
else
    
    save([output_dir filesep  output_name '_summary.mat'],'good_frames', 'all_subjects_minutes', 'all_mean_FD', 'all_poss_frames', 'all_possible_minutes', 'all_motion_conc_file')
end

disp('Done')
