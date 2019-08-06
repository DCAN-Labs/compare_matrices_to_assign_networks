
%load conc file of test dscalars.
%A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztemplate_match_intervals/MSC_sofar.conc');
%A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/infomap_all_subjects1run.conc');
%A = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/Community_Detection_Min_Dist_20_TieDen_0.02_MinNet_Size_400_MinReg_Size_30/infomap_1run.conc');%load confile of reference dscalars.
%B = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half2.conc');
%B= importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/infomap_all_subjects1run.conc');
template_matching = 1;
if template_matching == 1
%A = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/Zscored_dscalars/half1/intervals/MSChalf1_templ_matching_at_intervals.conc');
%B = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/Zscored_dscalars/all_frames_half2.conc');

A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/MSC_templ_match_intervals_dscalars.conc');
B = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half2.conc');

%E = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half1.conc');

else
end

infomap = 0;
if infomap == 1
C = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/MSC_infomap_intervals_dscalars.conc');

%C = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/Zscored_dscalars/intervals/MSChalf1_infomap_at_intervals.conc');
%D = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/MSC_Exacloud_lustre_backup/Infomap/Z-scoredconn_communities_half2/all_frames_half2.conc');
D = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/Infomap/Z-scoredconn_communities_half2/MSC_half2_allframes_infomap.conc');
%F = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/Infomap/Z-scoredconn_communities_half1/MSC_half1_allframes_infomap.conc');
else
end

if infomap == 0 && template_matching == 0
    disp('Neither infomap or template matching methods selected.  Exiting...')
    return
else
end

%subject should be arranged in following way in the conc file:
%subject1_time1_rep1
%subject1_time1_rep2
%subject1_time2_rep1
%subject1_time2_rep2
%subject2_time1_rep1
%etc

%parameters:
num_subjects = 10;
num_time_intervals = 8;
none_minuteslimit = 0;
num_reps = 10;
minutes = [1 2 3 4 5 10 15 20];

%% Adding paths for this function
this_code = which('plot_template_matching_reps');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

disp('Attempting to add neccesaary paths and functions.')
warning('off') %supress addpath warnings to nonfolders.
for i=2:np
    addpath(genpath(settings.path{i}));
end
%rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')

wb_command=settings.path_wb_c; %path to wb_command

%% Do checks
if template_matching == 1
    for i = 1:length(A)
        if exist(A{i}, 'file') == 0
            disp(['NOTE: Subject Series ' num2str(i) ' does not exist for template matching']);
            return
        else
        end
    end
    disp('All template matching test scalars files exist continuing ...')
    
    for i = 1:length(B)
        if exist(B{i}, 'file') == 0
            disp(['NOTE: Subject Series reference scalar ' num2str(i) ' does not exist for template matching']);
            return
        else
        end
    end
    disp('All template matching reference scalars files exist continuing ...')
    
else
end

if infomap == 1
    for i = 1:length(C)
        if exist(C{i}, 'file') == 0
            disp(['NOTE: Subject Series ' num2str(i) ' does not exist for infomap']);
            return
        else
        end
    end
    disp('All infomap test scalars files exist continuing ...')
    
    for i = 1:length(D)
        if exist(D{i}, 'file') == 0
            disp(['NOTE: Subject Series reference scalar  ' num2str(i) ' does not exist for infomap']);
            return
        else
        end
    end
    disp('All infomap reference scalars files exist continuing ...')
    
else
end

%% check minutes limit for all frames
if template_matching == 1
    if none_minuteslimit == 1
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects) + num_subjects;
    else
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects);
    end
    
    if size(A,1) == expected_num_scalars
        disp('Number of subjects and reps matches conc file')
    else
        disp('Number of subjects/reps is calculation does not match conc file.')
        return
    end
    
    
    if size(B,1) ~= num_subjects
        disp('Number of subjects in the conc file does not match the number of subjects')
        return
    else
        disp('Number of subjects in the conc file does matches the number of subjects....continuing')
    end
else
end

if infomap == 1
    
    if size(C,1) == expected_num_scalars
        disp('Number of subjects and reps matches conc file')
    else
        disp('Number of subjects/reps is calculation does not match conc file.')
        return
    end
    
    if size(D,1) ~= num_subjects
        disp('Number of subjects in the conc file does not match the number of subjects')
        return
    else
        disp('Number of subjects in the conc file does matches the number of subjects....continuing')
    end
else
end

if template_matching == 1
m = 1; % counter for line in conc file.
if none_minuteslimit == 0
    %% Step 4: load template networks from 2nd half of data:
    for i = 1:num_subjects
        reference_cii = ciftiopen(B{i},wb_command);
        other_half_networks = reference_cii.cdata;
        disp(['Reference scalar = ' B{i}])
        for j = 1:num_time_intervals
            for k = 1:num_reps
                test_cii = ciftiopen(A{m},wb_command);
                new_subject_labels = test_cii.cdata;
                disp(['Test scalar = ' A{m}]);
                disp('Caluating mutual information between scalars')
                %% Step 5: Calculate mutual information
                muI_templ(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                [VIn_templ(i,j,k), MIn_templ(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                m =m+1;
            end
        end
        disp(num2str(i))
    end
else
    %% Step 4: load template networks from 2nd half of data:
    for i = 1:num_subjects
        reference_cii = ciftiopen(B{i},wb_command);
        other_half_networks = reference_cii.cdata;
        disp(['Reference scalar = ' B{i}])
        for j = 1:num_time_intervals+1
            if j <= num_time_intervals %check if the interval is the last one (i.e. "none minutes limit" for which there will be no reps.)
                for k = 1:num_reps
                    test_cii = ciftiopen(A{m},wb_command);
                    new_subject_labels = test_cii.cdata;
                    disp(['Test scalar = ' A{m}]);
                    disp('Caluating mutual information between scalars')
                    %% Step 5: Calculate mutual information
                    muI_templ(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                    [VIn_templ(i,j,k), MIn_templ(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                    m =m+1;
                end
            else % if it is the minutes limit, dont go through reps.
                m =m+1;
            end
        end
        disp(num2str(i))
    end
end

else
end


if infomap == 1
m = 1; % counter for line in conc file.
if none_minuteslimit == 0
    %% Step 4: load template networks from 2nd half of data:
    for i = 1:num_subjects
        reference_cii = ciftiopen(D{i},wb_command);
        other_half_networks = reference_cii.cdata;
        disp(['Reference scalar = ' D{i}])
        for j = 1:num_time_intervals
            for k = 1:num_reps
                test_cii = ciftiopen(C{m},wb_command);
                new_subject_labels = test_cii.cdata;
                disp(['Test scalar = ' C{m}]);
                disp('Caluating mutual information between scalars')
                %% Step 5: Calculate mutual information
                muI_info(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                [VIn_info(i,j,k), MIn_info(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                m =m+1;
            end
        end
        disp(num2str(i))
    end
else
    %% Step 4: load template networks from 2nd half of data:
    for i = 1:num_subjects
        reference_cii = ciftiopen(D{i},wb_command);
        other_half_networks = reference_cii.cdata;
        disp(['Reference scalar = ' D{i}])
        for j = 1:num_time_intervals+1
            if j <= num_time_intervals %check if the interval is the last one (i.e. "none minutes limit" for which there will be no reps.)
                for k = 1:num_reps
                    test_cii = ciftiopen(C{m},wb_command);
                    new_subject_labels = test_cii.cdata;
                    disp(['Test scalar = ' C{m}]);
                    disp('Caluating mutual information between scalars')
                    %% Step 5: Calculate mutual information
                    muI_info(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                    [VIn_info(i,j,k), MIn_info(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                    m =m+1;
                end
            else % if it is the minutes limit, dont go through reps.
                m =m+1;
            end
        end
        disp(num2str(i))
    end
end

else
end
if infomap ==1
save('MSC_templ_info_10reps_MuI_to_2nd_half','muI_info','muI_templ','MIn_templ','MIn_info')
else
save('MSC_templ_10reps_MuI_to_2nd_half','muI_templ','MIn_templ')
end
%rep1 = squeeze(muI(:,:,1)');
%plot(minutes, rep1);
%C = permute(muI,[1,2,3]);
%D = reshape(C,[],size(muI,2),1);
%groups = [[1:10]'; [1:10]'; [1:10]';];
%groups = num2str([[1:10]'; [1:10]'; [1:10]';]);
%subjects = cellstr(groups);
%Dgroups = repmat([1:10]',1,8,3);
%rep1b = squeeze(muI(:,:,1)); rep2b = squeeze(muI(:,:,2)); rep3b = squeeze(muI(:,:,3));E = [rep1b;rep2b;rep3b;];
%parallelcoords(E,'group',subjects,'labels',timelabels,'quantile',0.25)
if template_matching == 1
mean_MuI_templ = mean(muI_templ,3); % calculate the mean of each interval (For each participant) using the reps.
SEM_MuI_templ = std(muI_templ,0,3)/sqrt(size(muI_templ,3)); % calculate the standard devitation of each interval (For each participant) using the reps.
else
end

if infomap == 1
mean_MuI_info = mean(muI_info,3); % calculate the mean of each interval (For each participant) using the reps.
SEM_MuI_info = std(muI_info,0,3)/sqrt(size(muI_info,3)); % calculate the standard devitation of each interval (For each participant) using the reps.
else
end

%% Actually make plots
if template_matching == 1
subplot (1,2,1)
for i = 1:num_subjects
errorbar(minutes,mean_MuI_templ(i,:),SEM_MuI_templ(i,:),'LineWidth',2); hold on
end
title('Template matching: Mutual information of various intervals to "hold out" half');
legend('Location','southeast');axis([0 20 0 2.4])
subplot (1,2,2)
else
end

if infomap == 1
for i = 1:num_subjects
errorbar(minutes,mean_MuI_info(i,:),SEM_MuI_templ(i,:),'LineWidth',2); hold on
end
title('Infomap: Mutual information of various intervals to "hold out" half');    
legend('Location','southeast');axis([0 20 0 2.4])    
else
end
disp('Done running plot_templatematching_reps');
set(gcf,'color','w')

