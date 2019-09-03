function check_twins_motion(all_motion_conc_file, TR, FD_column)

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


conc = strsplit(all_motion_conc_file, '.');
conc = char(conc(end));

if strcmp('conc',conc) == 1
    all_motion_conc = importdata(all_motion_conc_file);
else
    all_motion_conc = {all_motion_conc_file};
end

%The number 31 here means FD =0.3, for FD=0.2 use 21
%TR = 2.2;

for i = 1:length(all_motion_conc)
    [~,~,c] = fileparts(all_motion_conc{1});
    
    if strcmp(c,'.txt')
        FD_vector = importdata(all_motion_conc{i});
        good_frames(i) = sum(FD_vector);
    else
        
        
        load(all_motion_conc{i})
        FD02_vector =  motion_data{21}.frame_removal; % use 21 for FD= 0.2, 31 for 0.3
        good_frames(i) = sum(abs(FD02_vector-1));
        disp(i);
    end
end

disp('Number of good frames;')
good_frames';

disp('Amount of good data in minutes')



for i = 1:length(all_motion_conc)
    if strcmp(c,'.txt') % if TR is available from .mat file, use it to calculate the minutes of good data.
        good_frames_in_minutes=good_frames(i)*TR/60;
        poss_frames = length(importdata(all_motion_conc{i}));
        possible_minutes = poss_frames*TR/60;
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
    else
        load(all_motion_conc{i})
        good_frames_in_minutes=good_frames(i)*motion_data{1,FD_column}.epi_TR/60;
        poss_frames =  length(motion_data{FD_column}.frame_removal);
        possible_minutes = poss_frames*motion_data{1,FD_column}.epi_TR/60;
        disp([num2str(good_frames_in_minutes) ' minutes available out of ' num2str(possible_minutes) ' possible minutes'])
        
    end
end
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
        FD_vector =  motion_data{21}.frame_removal;
        for k = 1:splits %3 splits
            restsize = round(length(FD_vector)/splits);
            Rest(:,k)=FD_vector((restsize*k-restsize+1):(restsize*k));
            if exist('save_split_motion_vectors,var') == 1 && save_split_motion_vectors ==1
            else
                %save a list of files being written
                allrestlist{j} = [fileparts(fileparts(fileparts(lots_o_data_subs{i}))) '/MNINonLinear/Results/REST' num2str(k) '/FNL_preproc/REST_split_FD_mask.txt'];
                %save files
                %dlmwrite([fileparts(fileparts(fileparts(lots_o_data_subs{i}))) '/MNINonLinear/Results/REST' num2str(k) '/FNL_preproc/REST_split_FD_mask.txt'],Rest(:,k))
            end
        end
    end
else
end
disp('Done')