motion_all = importdata('/mnt/max/shared/projects/NIGGTWINS/Experiments/Template_matching/ADHD_Twins_12pairs_motion_file.conc');
for i = 1:length(motion_all)
    load(motion_all{i})
    FD02_vector =  motion_data{21}.frame_removal;
    good_frames(i) = sum(abs(FD02_vector-1));
end
good_frames