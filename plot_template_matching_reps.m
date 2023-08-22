
%% Description

%This code calculaters the normalized mutual information and a set od dscalars that were created with various time intervals, with an output dscalar.

%subject should be arranged in following way in the conc file:
%subject1_time1_rep1
%subject1_time1_rep2
%subject1_time2_rep1
%subject1_time2_rep2
%subject2_time1_rep1
%etc

%load conc file of test dscalars.
%A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztemplate_match_intervals/MSC_sofar.conc');
%A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/infomap_all_subjects1run.conc');
%A = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/Community_Detection_Min_Dist_20_TieDen_0.02_MinNet_Size_400_MinReg_Size_30/infomap_1run.conc');%load confile of reference dscalars.
%B = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half2.conc');
%B= importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/infomap_all_subjects1run.conc');

interval_setA = 1;
interval_setB = 1;

%-----------------------------------------------
if interval_setA == 1
    %A = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/Zscored_dscalars/half1/intervals/MSChalf1_templ_matching_at_intervals.conc');
    %B = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/Zscored_dscalars/all_frames_half2.conc');
    
    %A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/MSC_templ_match_intervals_dscalars.conc');
    %A = importdata('/home/faird/shared/projects/MSC_to_DCAN/split_halves/half1/MSC_half1_intervals_continuous.conc');
    %B = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half2.conc');
    A = importdata('/home/rando149/hermosir/subpop_subjects/conc_files/template_matching_reli_surface_only.conc');
    %B = importdata('/home/faird/shared/projects/MSC_to_DCAN/analyses/template_matching/all_frames/MSC_to_DCAN_all_frames_template_matching_half2.conc');
    %B = importdata('/home/faird/shared/projects/MSC_to_DCAN/split_halves/half2ADHD315/MSC_half2_to_ADHD315Z_Zscored_recolored_dscalars.conc');
    %E = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half1.conc');
    B = importdata('/home/rando149/hermosir/subpop_subjects/conc_files/template_matching_nets_surface_only_halfb.conc');
else
end

if interval_setB == 1
    %C = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/MSC_infomap_intervals_dscalars.conc');
    
    %C = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/Zscored_dscalars/intervals/MSChalf1_infomap_at_intervals.conc');
    %D = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/MSC_Exacloud_lustre_backup/Infomap/Z-scoredconn_communities_half2/all_frames_half2.conc');
    %D = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/Infomap/Z-scoredconn_communities_half2/MSC_half2_allframes_infomap.conc');
    %F = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/Infomap/Z-scoredconn_communities_half1/MSC_half1_allframes_infomap.conc');
    C = importdata('/home/faird/shared/projects/MSC_to_DCAN/analyses/template_matching/full_intervals/MSChalf1_templ_matching_at_intervals_MSI_paths.conc');
    %D = importdata('/home/faird/shared/projects/MSC_to_DCAN/analyses/template_matching/all_frames/MSC_to_DCAN_all_frames_template_matching_half2.conc');
    D = importdata('/home/faird/shared/projects/MSC_to_DCAN/split_halves/half2ADHD315/MSC_half2_to_ADHD315Z_Zscored_recolored_dscalars.conc');
    
else
end

if interval_setB == 0 && interval_setA == 0
    disp('Neither infomap or template matching methods selected.  Exiting...')
    return
else
end


%parameters:
%num_subjects = 5;
num_subjects_setA = 5;
num_subjects_setB = 10;

use_cortex_only =1; %will be applied to both datasets.
num_time_intervals = 8;
none_minuteslimit = 0;
num_reps = 10;
minutes = [1 2 3 4 5 10 15 20];

subject_labelsA ={'sub01','sub02','sub03','sub04','sub05'};
subject_labelsB ={'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};

%% Adding paths for this function
this_code = which('plot_template_matching_reps');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

disp('Attempting to add neccesaary paths and functions.')
%warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end
%rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
%warning('on')

wb_command=settings.path_wb_c; %path to wb_command
reset(RandStream.getGlobalStream,sum(100*clock)); % added to ensure randomness of frame sampling during parellization -RH 05/30/2023

%% Do checks
if interval_setA == 1
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

if interval_setB == 1
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
if interval_setA == 1
    if none_minuteslimit == 1
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects_setA) + num_subjects_setA;
    else
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects_setA);
    end
    
    if size(A,1) == expected_num_scalars
        disp('Number of subjects and reps matches conc file')
    else
        disp('Number of subjects/reps is calculation does not match conc file.')
        return
    end
    
    
    if size(B,1) ~= num_subjects_setA
        disp('Number of subjects in the conc file does not match the number of subjects')
        return
    else
        disp('Number of subjects in the conc file does matches the number of subjects....continuing')
    end
else
end

if interval_setB == 1
    if none_minuteslimit == 1
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects_setB) + num_subjects_setB;
    else
        expected_num_scalars = (num_time_intervals * num_reps * num_subjects_setB);
    end
    
    if size(C,1) == expected_num_scalars
        disp('Number of subjects and reps matches conc file')
    else
        disp('Number of subjects/reps is calculation does not match conc file.')
        return
    end
    
    if size(D,1) ~= num_subjects_setB
        disp('Number of subjects in the conc file does not match the number of subjects')
        return
    else
        disp('Number of subjects in the conc file does matches the number of subjects....continuing')
    end
else
end

if interval_setA == 1
    m = 1; % counter for line in conc file.
    if none_minuteslimit == 0
        %% Step 4: load template networks from 2nd half of data:
        for i = 1:num_subjects_setA
            reference_cii = ciftiopen(B{i},wb_command);
            other_half_networks = reference_cii.cdata;
            if use_cortex_only ==1
                other_half_networks=other_half_networks(1:59412,1); % hardcode warning.
            end
            disp(['Reference scalar = ' B{i}])
            for j = 1:num_time_intervals
                for k = 1:num_reps
                    test_cii = ciftiopen(A{m},wb_command);
                    
                    new_subject_labels = test_cii.cdata;
                    if use_cortex_only ==1
                        new_subject_labels=new_subject_labels(1:59412,1); % hardcode warning.
                    end
                    disp(['Test scalar = ' A{m}]);
                    disp('Caluating mutual information between scalars')
                    %% Step 5: Calculate mutual information
                    muI_setA(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                    [VIn_setA(i,j,k), MIn_setA(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                    m =m+1;
                end
            end
            disp(num2str(i))
        end
    else
        %% Step 4: load template networks from 2nd half of data:
        for i = 1:num_subjects_setA
            reference_cii = ciftiopen(B{i},wb_command);
            other_half_networks = reference_cii.cdata;
            if use_cortex_only ==1
                other_half_networks=other_half_networks(1:59412,1); % hardcode warning.
            end
            
            disp(['Reference scalar = ' B{i}])
            for j = 1:num_time_intervals+1
                if j <= num_time_intervals %check if the interval is the last one (i.e. "none minutes limit" for which there will be no reps.)
                    for k = 1:num_reps
                        test_cii = ciftiopen(A{m},wb_command);
                        new_subject_labels = test_cii.cdata;
                        if use_cortex_only ==1
                            new_subject_labels=new_subject_labels(1:59412,1); % hardcode warning.
                        end
                        
                        disp(['Test scalar = ' A{m}]);
                        disp('Caluating mutual information between scalars')
                        %% Step 5: Calculate mutual information
                        muI_setA(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                        [VIn_setA(i,j,k), MIn_setA(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
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


if interval_setB == 1
    m = 1; % counter for line in conc file.
    if none_minuteslimit == 0
        %% Step 4: load template networks from 2nd half of data:
        for i = 1:num_subjects_setB
            reference_cii = ciftiopen(D{i},wb_command);
            other_half_networks = reference_cii.cdata;
            if use_cortex_only ==1
                other_half_networks=other_half_networks(1:59412,1); % hardcode warning.
            end
            disp(['Reference scalar = ' D{i}])
            for j = 1:num_time_intervals
                for k = 1:num_reps
                    test_cii = ciftiopen(C{m},wb_command);
                    
                    new_subject_labels = test_cii.cdata;
                    
                    if use_cortex_only ==1
                        new_subject_labels=new_subject_labels(1:59412,1); % hardcode warning.
                    end
                    disp(['Test scalar = ' C{m}]);
                    disp('Caluating mutual information between scalars')
                    %% Step 5: Calculate mutual information
                    muI_setB(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                    [VIn_setB(i,j,k), MIn_setB(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
                    m =m+1;
                end
            end
            disp(num2str(i))
        end
    else
        %% Step 4: load template networks from 2nd half of data:
        for i = 1:num_subjects_setB
            reference_cii = ciftiopen(D{i},wb_command);
            other_half_networks = reference_cii.cdata;
            if use_cortex_only ==1
                other_half_networks=other_half_networks(1:59412,1); % hardcode warning.
            end
            disp(['Reference scalar = ' D{i}])
            for j = 1:num_time_intervals+1
                if j <= num_time_intervals %check if the interval is the last one (i.e. "none minutes limit" for which there will be no reps.)
                    for k = 1:num_reps
                        test_cii = ciftiopen(C{m},wb_command);
                        
                        new_subject_labels = test_cii.cdata;
                        if use_cortex_only ==1
                            new_subject_labels=new_subject_labels(1:59412,1); % hardcode warning.
                        end
                        disp(['Test scalar = ' C{m}]);
                        disp('Caluating mutual information between scalars')
                        %% Step 5: Calculate mutual information
                        muI_setB(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                        [VIn_setB(i,j,k), MIn_info(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
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

% if infomap ==1
%     save('MSC_templ_info_10reps_MuI_to_2nd_half','muI_info','muI_templ','MIn_templ','MIn_info')
% else
%     save('MSC_templ_10reps_MuI_to_2nd_half','muI_templ','MIn_templ')
% end

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
if interval_setA == 1
    %mean_MuI_templ = mean(muI_setA,3); % calculate the mean of each interval (For each participant) using the reps.
    %SEM_MuI_templ = std(muI_setA,0,3)/sqrt(size(muI_setA,3)); % calculate the standard devitation of each interval (For each participant) using the reps.
    mean_MIn_setA = mean(MIn_setA,3);
    SEM_MIn_setA = std(MIn_setA,0,3)/sqrt(size(MIn_setA,3));
    
else
end

if interval_setB == 1
    %mean_MuI_info = mean(muI_info,3); % calculate the mean of each interval (For each participant) using the reps.
    %SEM_MuI_info = std(muI_info,0,3)/sqrt(size(muI_info,3)); % calculate the standard devitation of each interval (For each participant) using the reps.
    mean_MIn_setB = mean(MIn_setB,3);
    SEM_MIn_setB = std(MIn_setB, 0,3)/sqrt(size(MIn_setB,3));
    
else
end

%% Actually make plots

% if template_matching == 1
%     subplot (1,2,1)
%     for i = 1:num_subjects
%         errorbar(minutes,mean_MuI_templ(i,:),SEM_MuI_templ(i,:),'LineWidth',2); hold on
%     end
%     title('Template matching: Mutual information of various intervals to "hold out" half');
%     legend('Location','southeast');axis([0 20 0 2.4])
%     subplot (1,2,2)
% else
% end
%
% if infomap == 1
%     for i = 1:num_subjects
%         errorbar(minutes,mean_MuI_info(i,:),SEM_MuI_templ(i,:),'LineWidth',2); hold on
%     end
%     title('Infomap: Mutual information of various intervals to "hold out" half');
%     legend('Location','southeast');axis([0 20 0 2.4])
% else
% end



F = figure();
set(gcf,'color','w')
set(gca,'FontSize',16)

if interval_setA == 1
    
    if interval_setA ==1 && interval_setB ==1
        subplot (1,2,1)
    end
    %HARDCODE warning
    % MSC
    %     max_null_MIn = 0.4004;
    %     average_null_MIn = 0.3663;
    %     min_null_MIn =0.3208;
    
    
    %Subpop
    max_null_MIn = 0.5813;
    average_null_MIn = 0.4764;
    min_null_MIn =0.4343;
    
    patch([0 20 20 0],[min_null_MIn min_null_MIn max_null_MIn max_null_MIn],[180/255 180/255 180/255],'FaceAlpha',0.3); hold on
    line([0 20],[max_null_MIn,max_null_MIn],'LineStyle','--');hold on
    line([0 20],[average_null_MIn,average_null_MIn],'Color','black','LineWidth',3,'LineStyle',':');hold on
    line([0 20],[min_null_MIn,min_null_MIn],'LineStyle','--');hold on
    
    %MSC
    %     max_real_MIn = 0.6891;
    %     average_real_MIn = 0.6101;
    %     min_real_MIn = 0.5136;
    %subpop
    max_real_MIn = 0.7963;
    average_real_MIn = 0.7571;
    min_real_MIn = 0.6995;
    
    patch([0 20 20 0],[min_real_MIn min_real_MIn max_real_MIn max_real_MIn],[0/255 120/255 200/255],'FaceAlpha',0.3); hold on
    line([0 20],[max_real_MIn,max_real_MIn],'LineStyle','--'); hold on
    line([0 20],[average_real_MIn,average_real_MIn],'Color','blue','LineWidth',3,'LineStyle',':');hold on
    line([0 20],[min_real_MIn,min_real_MIn],'LineStyle','--');hold on
    
    for i = 1:num_subjects_setA
        G(i) = errorbar(minutes,mean_MIn_setA(i,:),SEM_MIn_setA(i,:),'LineWidth',2); hold on
    end
    
    %MSC_labels ={'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};
    %lgd = legend(G,MSC_labels,'Location','southeast');
    lgd = legend(G,subject_labelsA,'Location','southeast');
    lgd.FontSize = 10;
    %axis([0 20 0.25 0.70]);
    %axis([0 20 0.1 0.70]);
    axis([0 20 0.1 0.8]);
end

if interval_setB == 1
    
    if interval_setA ==1 && interval_setB ==1
        subplot (1,2,2)
    end
    
    %HARDCODE warning
    % MSC
    %     max_null_MIn = 0.4004;
    %     average_null_MIn = 0.3663;
    %     min_null_MIn =0.3208;
  
    % MSC
         max_null_MIn = 0.4629;
         average_null_MIn = 0.4247;
         min_null_MIn =0.3755;   
    
    %Subpop
    %max_null_MIn = 0.4629;
    %average_null_MIn = 0.4247;
    %min_null_MIn =0.4343;
    
    patch([0 20 20 0],[min_null_MIn min_null_MIn max_null_MIn max_null_MIn],[180/255 180/255 180/255],'FaceAlpha',0.3); hold on
    line([0 20],[max_null_MIn,max_null_MIn],'LineStyle','--');hold on
    line([0 20],[average_null_MIn,average_null_MIn],'Color','black','LineWidth',3,'LineStyle',':');hold on
    line([0 20],[min_null_MIn,min_null_MIn],'LineStyle','--');hold on
    
    %MSC w cortex
    %      max_real_MIn = 0.6891;
    %      average_real_MIn = 0.6101;
    %      min_real_MIn = 0.5136;
    
    %subpop
    %     max_real_MIn = 0.7963;
    %     average_real_MIn = 0.7571;
    %     min_real_MIn = 0.6995;
    %MSC_ cortex only
    max_real_MIn = 0.7356;
    average_real_MIn = 0.6790;
    min_real_MIn = 0.6072;
    
    patch([0 20 20 0],[min_real_MIn min_real_MIn max_real_MIn max_real_MIn],[0/255 120/255 200/255],'FaceAlpha',0.3); hold on
    line([0 20],[max_real_MIn,max_real_MIn],'LineStyle','--'); hold on
    line([0 20],[average_real_MIn,average_real_MIn],'Color','blue','LineWidth',3,'LineStyle',':');hold on
    line([0 20],[min_real_MIn,min_real_MIn],'LineStyle','--');hold on
    
    for i = 1:num_subjects_setB
        G(i) = errorbar(minutes,mean_MIn_setB(i,:),SEM_MIn_setB(i,:),'LineWidth',2); hold on
    end
    
    %MSC_labels ={'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};
    lgd = legend(G,subject_labelsB,'Location','southeast');
    lgd.FontSize = 10;
    %axis([0 20 0.25 0.70]);
    %axis([0 20 0.1 0.70]);
    axis([0 20 0.1 0.8]);
end

set(gcf,'color','w')
print(['NMI_intervals_' num2str(randi([1 1000],1))],'-dpng','-r300'); % provide a random number so that you don't overwrite figures. good idea? -RH

disp('Done running plot_templatematching_reps');

