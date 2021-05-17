%load('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/surface_area/lucianov_test_cleaned.mat');%load surface and volume data
%load('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/surface_area/ADHDsymptoms.mat')

%niggtwins
%dscalarwithassignments = importdata('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/template_matching_dscalars/template_matching_cleaned_dscalar.conc');
%HCP twins
%dscalarwithassignments = importdata('/mnt/max/shared/projects/NIGGTWINS/WTO/Data/HCP_data/monozygotic_twins_100pairs.conc');

%ABCD twins
%dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_32_rtg_monozyg.conc');
%dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_180_rtg_dizyg.conc');
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_twins_surface_area/ABCD_38_rtg_monozyg_ID_10min.conc');

%MSC Halves
%templ_match_dscalarwithassignments = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/bothhalvesdsclars.conc');
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/template_matching/all_frames/MSC_to_DCAN_all_frames_templ.conc');
%infomap_dscalarwithassignments = importdata('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/MSC_Exacloud_lustre_backup/Infomap/bothhalvesdsclars.conc');
%infomap_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/infomap_scott/infomap_dscalars.conc');

%cross-method comparison for infomap and template matching
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/cross-method/MSC_to_DCAN_all_frames_templ_match_half1_infomap_half1.conc');
%infomap_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/cross-method/MSC_to_DCAN_all_frames_templ_match_half2_infomap_half2.conc');

%within half analysis
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/cross-method/MSC_to_DCAN_all_frames_templ_match_half1_infomap_half2.conc');
%infomap_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/cross-method/MSC_to_DCAN_all_frames_templ_match_half2_infomap_half1.conc');

%overlapping MSC
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/template_matching/all_frames/MSC_to_DCAN_all_frames_overlapping_templ.conc');
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/MSC_to_DCAN/analyses/template_matching/all_frames/MSC_to_DCAN_all_frames_overlapping_templ_recolored.conc');
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/best10_ABCDsubs/split_halves/best10_ABCDsubs_TM_overlap_dscalar.conc');
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/best10_ABCDsubs/split_halves/best10_ABCDsubs_TM_overlap_dtseries.conc');

%ABCD_GROUP_AVG
%templ_match_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/ABCD_GROUP_AVERAGES/templ_only.conc');

%ABCD_best_subs
templ_match_dscalarwithassignments = importdata('/home/faird/shared/projects/ABCD_net_template_matching/best10_ABCDsubs/split_halves/best10_ABCDsubs_TM_singlenet_dscalar.conc');
%infomap_dscalarwithassignments = importdata('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/best10_ABCDsubs/split_halves/infomap/merged_densities/infomap_both_halves.conc');

twins=0;
surface_only =0; ncortgrey = 59412;
omnibus =0; %for overlapping networks, set to true if you want to concatenate all networks together and run an omnibus mutual information.  Otherwise MuI will be calculated for each network seperately.
%note: network names listed below have empty networks removed (e.g. 4 and 6 are revmoed from the data.)
network_names = {   'DMN'    'Vis'    'FP'      'DAN'       'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};

%% Add necessary paths
addpath ('/panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks')
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools'))

this_code = which('pairwise_Mutualinfo');
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
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')

wb_command=settings.path_wb_c; %path to wb_command
tic

%check to make sure that surface files exist
if exist('templ_match_dscalarwithassignments','var') == 1
    for i = 1:length(templ_match_dscalarwithassignments)
        if rem(i,100)==0
            disp([' Validating file existence ' num2str(i)]);toc;
        end
        if exist(templ_match_dscalarwithassignments{i},'file') == 0
            NOTE = ['Subject dscalar ' num2str(i) ' does not exist']
            disp(templ_match_dscalarwithassignments{i});
            return
        else
        end
    end
    disp('All template matching dscalars exist continuing ...')
end

if exist('infomap_dscalarwithassignments','var') == 1
    for i = 1:length(infomap_dscalarwithassignments)
        if rem(i,100)==0
            disp([' Validating file existence ' num2str(i)]);toc;
        end
        if exist(infomap_dscalarwithassignments{i},'file') == 0
            NOTE = ['Subject dscalar ' num2str(i) ' does not exist']
            disp(infomap_dscalarwithassignments{i});
            return
        else
        end
    end
    disp('All infomap dscalars exist continuing ...')
end

if twins == 1
    Number_subjects = length(dscalarwithassignments);
    network_assignment_filetype = strsplit(dscalarwithassignments{1}, '.');
else
    Number_subjects = length(templ_match_dscalarwithassignments);
    network_assignment_filetype = strsplit(templ_match_dscalarwithassignments{1}, '.');
end

cifti_type = char(network_assignment_filetype(end-1));
if strcmp('dtseries',cifti_type) == 1
    overlap =1;
else
    overlap =0;
end

disp('loading scalars')
if twins == 1
    for i = 1:Number_subjects
        allscalars_templ_twins = ciftiopen(dscalarwithassignments{i},wb_command);
        if overlap ==1
            thissubjectsgreys = single(allscalars_templ_temp.cdata);
            thissubjectsgreys(:,4) = []; % remove networks with no assingments
            thissubjectsgreys(:,5) = []; % remove networks with no assingments
            thissubjectsgreys  = thissubjectsgreys+1;
            if omnibus ==1
                thissubjectsgreys  = reshape(thissubjectsgreys,[],1); % reshape
            end
            %add 1 to all values for VIN calculation.  Does not like "0" %values.
            
            allscalars_templ(:,i) = single(thissubjectsgreys);
            
        else
            if exist('surface_only','var') == 1 && surface_only ==1
                thissubjectsgreys = single(allscalars_templ_twins.cdata);
                thissubjectsgreys = thissubjectsgreys(1:ncortgrey,:);
                allscalars_templ(:,i) = single(thissubjectsgreys);
            else
                allscalars_templ(:,i) = single(allscalars_templ_twins.cdata);
            end
        end
    end
    
else
    Number_subjects = length(templ_match_dscalarwithassignments);
    
    for i = 1:Number_subjects
        allscalars_templ_temp = ciftiopen(templ_match_dscalarwithassignments{i},wb_command);
        if overlap ==1
            thissubjectsgreys = single(allscalars_templ_temp.cdata);
            thissubjectsgreys(:,4) = []; % remove networks with no assingments
            thissubjectsgreys(:,5) = []; % remove networks with no assingments
            thissubjectsgreys  = thissubjectsgreys+1;
            if omnibus ==1
                thissubjectsgreys  = reshape(thissubjectsgreys,[],1); % reshape
            end
            %add 1 to all values for VIN calculation.  Does not like "0" %values.
            
            allscalars_templ(:,:,i) = single(thissubjectsgreys);
        else
            if exist('surface_only','var') == 1 && surface_only ==1
                thissubjectsgreys = single(allscalars_templ_temp.cdata);
                thissubjectsgreys = thissubjectsgreys(1:ncortgrey,:);
                allscalars_templ(:,i) = single(thissubjectsgreys);
            else
                allscalars_templ(:,i) = single(allscalars_templ_temp.cdata);
            end
        end
    end
    
    if exist('infomap_dscalarwithassignments','var') == 1
        for i = 1:Number_subjects
            allscalars_info_temp = ciftiopen(infomap_dscalarwithassignments{i},wb_command);
            if overlap ==1
                thissubjectsgreys = single(allscalars_info_temp.cdata);
                thissubjectsgreys(:,4) = []; % remove networks with no assingments
                thissubjectsgreys(:,5) = []; % remove networks with no assingments
                thissubjectsgreys  = thissubjectsgreys+1;
                if omnibus ==1
                    thissubjectsgreys  = reshape(thissubjectsgreys,[],1); % reshape
                end
                %add 1 to all values for VIN calculation.  Does not like "0" %values.
                
                allscalars_info(:,i) = single(thissubjectsgreys);
            else
                if exist('surface_only','var') == 1 && surface_only ==1
                    thissubjectsgreys = single(allscalars_info_temp.cdata);
                    thissubjectsgreys = thissubjectsgreys(1:ncortgrey,:);
                    allscalars_info(:,i) = single(thissubjectsgreys);
                else
                    allscalars_info(:,i) = single(allscalars_info_temp.cdata);
                end
            end
        end
    else
    end
end

if twins == 1
    twin1_indices = 1:2:length(dscalarwithassignments);
    twin2_indices = 2:2:length(dscalarwithassignments);
else
    twin1_indices = 1:2:length(templ_match_dscalarwithassignments);
    twin2_indices = 2:2:length(templ_match_dscalarwithassignments);
end


%  for j = 1:16 %go through every networks
%      for i = 1:26 %go through every subject
%          for k = 1:26 %go through every other subject
%              if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
%                  %if i < k
%                 subject1 = network_surfarea (i,j); eta_net_assign1 = allscalars(i);
%                  subjectk = network_surfarea (k,j); eta_net_assignk = allscalars(k);
%                  all_network_difference_surfarea(i,k,j) = abs(subject1 - subjectk);
%              end
%          end
%      end
%  end


% for i=1:16   %go through every network
% diffmatrix = all_network_difference_surfarea(:,:,i);
% diff_indices = find(diffmatrix~=0);
% diff_values = diffmatrix(diff_indices);
% a = mean(network_difference_surfarea(:,i));
% proportion_of_diff_below_paired_average_diff(i) = length(find(diffmatrix(diff_indices)<a))/length(diff_indices);
% %histogram(diff_values,30); hold on; histogram(network_difference_surfarea(:,i),30,'FaceColor','r');
% %hist(diffmatrix(diff_indices),30); hold on; bar(a,16,50,'r')
% end
if overlap ==1
    V = size(allscalars_templ,2);
else
    V = 1;
end

for j=1:V
    disp(['Network ' num2str(j)])
    %calculate mutual information for real pairs
    disp('Calculating mutual information for real pairs')
    for i = 1:2:Number_subjects %number of subjects
        if overlap ==1
            muI_templ(round(i/2),j) = MutualInformation(allscalars_templ(:,j,i),allscalars_templ(:,j,i+1)); %Mutual information
            [VIn_templ(round(i/2),j), MIn_templ(round(i/2),j)] = partition_distance(allscalars_templ(:,j,i),allscalars_templ(:,j,i+1)); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
        else
            muI_templ(round(i/2),j) = MutualInformation(allscalars_templ(:,i),allscalars_templ(:,i+1)); %Mutual information
            [VIn_templ(round(i/2),j), MIn_templ(round(i/2),j)] = partition_distance(allscalars_templ(:,i),allscalars_templ(:,i+1)); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            
        end
        
    end
    
    
    %calculate mututal information to all possible combinations, except pairs
    disp('Calculating mutual information to all possible combinations, except pairs')
    for i = 1:Number_subjects %go through every subject
        for k = 1:Number_subjects %go through every other subject
            if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
                if overlap ==1
                    
                    eta_net_assign1 = allscalars_templ(:,j,i);
                    eta_net_assignk = allscalars_templ(:,j,k);
                else
                    eta_net_assign1 = allscalars_templ(:,i);
                    eta_net_assignk = allscalars_templ(:,k);
                end
                unpaired_muI(i,k,j) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                [unpaired_VIn(i,k,j), unpaired_MIn(i,k,j)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            end
        end
    end
    
    % Calculate mututal information to from twin1 to all twin2s.  Just used for
    % graph.
    m = 1;n=1;
    disp('Calculating mutual information for twin groups for graph')
    for i = twin1_indices %go through every subject
        for k = twin2_indices %go through every other subject
            %if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
            if overlap ==1
                
                eta_net_assign1 = allscalars_templ(:,j,i);
                eta_net_assignk = allscalars_templ(:,j,k);
            else
                eta_net_assign1 = allscalars_templ(:,i);
                eta_net_assignk = allscalars_templ(:,k);
            end
            block_muI_templ(m,n,j) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
            [block_VIn_templ(m,n,j), block_MIn_templ(m,n,j)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            %end
            n=n+1;
        end
        n=1;m=m+1;
    end
    
    
    %calculate mututal informatin to all possible combinations, including pairs
    disp('Calculating mutual information to all possible combinations, including pairs')
    for i = 1:Number_subjects %go through every subject
        for k = 1:Number_subjects %go through every other subject
            %if i < k% && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
            if overlap ==1
                
                eta_net_assign1 = allscalars_templ(:,j,i);
                eta_net_assignk = allscalars_templ(:,j,k);
            else
                eta_net_assign1 = allscalars_templ(:,i);
                eta_net_assignk = allscalars_templ(:,k);
            end
            allposs_muI_pairs(i,k,j) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
            [allposs_w_pairs_VIn(i,k,j), allposs_w_pairs_MIn(i,k,j)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
            %end
        end
    end
    
    
    diff_muI_matrix(:,:,j) = unpaired_muI(:,:,j);
    diff_muI_indices = find(unpaired_muI(:,:,j)~=0);
    this_net_muI_matrix = diff_muI_matrix(:,:,j);
    diff_muI_values(:,j) = this_net_muI_matrix(diff_muI_indices);
    
    diff_VIn_matrix(:,:,j) = unpaired_VIn(:,:,j);
    diff_VIn_indices = find(unpaired_VIn(:,:,j)~=0);
    this_net_VIn_matrix = diff_VIn_matrix(:,:,j);
    diff_VIn_values(:,j) = this_net_VIn_matrix(diff_VIn_indices);
    
    diff_MIn_matrix(:,:,j) = unpaired_MIn(:,:,j);
    diff_MIn_indices = find(unpaired_MIn(:,:,j)~=0);
    this_net_MIn_matrix = diff_MIn_matrix(:,:,j);
    diff_MIn_values(:,j) = this_net_MIn_matrix(diff_MIn_indices);
    
    %t = (a - b) / (c/sqrt(length(diff_muI_indices)));
    %hist(diff_muI_values,30); hold on; bar(a,16,'BarWidth',0.01,'FaceColor','r');
    
    %figure()
    %plot(muI)
    %plot(muI,'MarkerSize',10,'MarkerEdgeColor','b')
    %boxplot(info_muI')
    if exist('infomap_dscalarwithassignments','var') == 1
        for i = 1:2:Number_subjects %number of subjects
            muI_info(round(i/2),1) = MutualInformation(allscalars_info(:,i),allscalars_info(:,i+1)); %Mutual information
            [VIn_info(round(i/2),1), MIn_info(round(i/2),1)] = partition_distance(allscalars_info(:,i),allscalars_info(:,i+1)); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
        end
        
        disp('Calculating mutual information to all possible combinations, except pairs')
        
        for i = 1:Number_subjects %go through every subject
            for k = 1:Number_subjects %go through every other subject
                if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
                    eta_net_assign1 = allscalars_info(:,i);
                    eta_net_assignk = allscalars_info(:,k);
                    unpaired_muI_info(i,k) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                    [unpaired_VIn_info(i,k), unpaired_MIn_info(i,k)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
                end
            end
        end
        
        m = 1;n=1;
        disp('Calculating mutual information for twin groups for graph')
        for i = twin1_indices %go through every subject
            for k = twin2_indices %go through every other subject
                %if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
                eta_net_assign1 = allscalars_info(:,i);
                eta_net_assignk = allscalars_info(:,k);
                block_muI_info(m,n) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                [block_VIn_info(m,n), block_MIn_info(m,n)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
                %end
                n=n+1;
            end
            n=1;m=m+1;
        end
        
        
        disp('Calculating mutual information to all possible combinations, including pairs')
        for i = 1:Number_subjects %go through every subject
            for k = 1:Number_subjects %go through every other subject
                %if i < k% && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
                eta_net_assign1 = allscalars_info(:,i);
                eta_net_assignk = allscalars_info(:,k);
                allposs_muI_pairs_info(i,k) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                [allposs_w_pairs_VIn_info(i,k), allposs_w_pairs_MIn_info(i,k)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
                %end
            end
        end
        
        diff_muI_matrix_info = unpaired_muI_info;
        %[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
        diff_muI_indices_info = find(unpaired_muI_info~=0);
        diff_muI_values_info = diff_muI_matrix_info(diff_muI_indices_info);
        %a = mean(abs(muI));b = mean(diff_muI_values); c = std(diff_muI_values);
        %proportion_of_mu_I_diff_below_paired_average_diff = length(find(diff_muI_values<a))/length(diff_muI_indices);
        
        diff_VIn_matrix_info = unpaired_VIn_info;
        %[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
        diff_VIn_indices_info = find(unpaired_VIn_info~=0);
        diff_VIn_values_info = diff_VIn_matrix_info(diff_VIn_indices_info);
        %a = mean(abs(VIn));b = mean(diff_VIn_values); c = std(diff_VIn_values);
        %proportion_of_mu_I_diff_below_paired_average_diff = length(find(diff_VIn_values<a))/length(diff_VIn_indices);
        
        diff_MIn_matrix_info = unpaired_MIn_info;
        %[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
        diff_MIn_indices_info = find(unpaired_MIn_info~=0);
        diff_MIn_values_info = diff_MIn_matrix_info(diff_MIn_indices_info);
        %a = mean(abs(MIn));b = mean(diff_MIn_values); c = std(diff_MIn_values);
        %proportion_of_mu_I_diff_below_paired_average_diff =
        %length(find(diff_MIn_values<a))/length(diff_MIn_indices);
    else
    end
end



figure()
if twins ==1
    
    subplot(2,4,1)
    imagesc(allposs_w_pairs_MIn);title('Template matching NMI matrix');
    subplot(2,4,2)
    %figure();
    histogram(diff_muI_values,30); hold on; histogram(muI_templ,30,'FaceColor','r');xlabel('Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real disc. twin pairs'},'Location','northeast'); title('Mutual info Shuffled pairs');
    subplot(2,4,3)
    %figure();
    histogram(diff_VIn_values,30); hold on; histogram(VIn_templ,30,'FaceColor','r');xlabel('Variation of Information (bits)');ylabel('count');legend({'All possible half pairs','Real disc. twin pairs'},'Location','northeast'); title('Variation of info Shuffled pairs');
    subplot(2,4,4)
    %figure();
    histogram(diff_MIn_values,30); hold on; histogram(MIn_templ,30,'FaceColor','r');xlabel('Normalized Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real disc. twin pairs'},'Location','northeast'); title('Norm. Mutual info Shuffled pairs');
    
else
end


subplot(2,4,1)
if overlap ==0
    imagesc(allposs_w_pairs_MIn);title('Template matching NMI matrix');
    subplot(2,4,2)
    %figure();
    histogram(diff_muI_values,30); hold on; histogram(muI_templ,30,'FaceColor','r');xlabel('Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Mutual info Shuffled pairs');
    subplot(2,4,3)
    %figure();
    histogram(diff_VIn_values,30); hold on; histogram(VIn_templ,30,'FaceColor','r');xlabel('Variation of Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Variation of info Shuffled pairs');
    subplot(2,4,4)
    %figure();
    histogram(diff_MIn_values,30); hold on; histogram(MIn_templ,30,'FaceColor','r');xlabel('Normalized Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Norm. Mutual info Shuffled pairs');
    
    if exist('infomap_dscalarwithassignments','var') == 1
        subplot(2,4,5)
        imagesc(allposs_w_pairs_MIn_info);title('Infomap NMI matrix');
        subplot(2,4,6)
        %figure();
        histogram(diff_muI_values_info,30); hold on; histogram(muI_info,30,'FaceColor','r');xlabel('Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Mutual info Shuffled pairs');
        subplot(2,4,7)
        %figure();
        histogram(diff_VIn_values_info,30); hold on; histogram(VIn_info,30,'FaceColor','r');xlabel('Variation of Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Variation of info Shuffled pairs');
        subplot(2,4,8)
        %figure();
        histogram(diff_MIn_values_info,30); hold on; histogram(MIn_info,30,'FaceColor','r');xlabel('Normalized Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast'); title('Norm. Mutual info Shuffled pairs');
    else
    end
else
end

figure()
if overlap ==1
    for b=1: size(allscalars_templ,2)
        subplot(4,4,b)
        this_block = squeeze(block_MIn_templ(:,:,b));
        imagesc(this_block);
    end
else
    imagesc(block_MIn_templ)
end

if exist('infomap_dscalarwithassignments','var') == 1
    figure()
    imagesc(block_MIn_info)
end

if twins == 1
    title('Normalized Mutual information from twin 1 to twin 2')
    xlabel('Twin2')
    ylabel('Twins1')
else
    title('Normalized Mutual information from half 1 to half 2')
    xlabel('Subject')
    ylabel('Subject')
end



% figure()
% subplot(1,2,1)
% G = graph(allposs_muI_pairs,{'MSC1a','MSC1b','MSC2a','MSC2b','MSC3a','MSC3b','MSC4a','MSC4b','MSC5a','MSC5b','MSC6a','MSC6b','MSC7a','MSC7b','MSC8a','MSC8b','MSC9a','MSC9b','MS10a','MSC10b'},'upper','omitselfloops');
% LWidths = abs((zscore(G.Edges.Weight/mean(G.Edges.Weight))));
% plot(G,'LineWidth',LWidths);
% disp('Done running pairwise_mutalinfo.')
% title('Zscored edge weights')
% subplot(1,2,2)
% H = graph(allposs_muI_pairs_info,{'MSC1a','MSC1b','MSC2a','MSC2b','MSC3a','MSC3b','MSC4a','MSC4b','MSC5a','MSC5b','MSC6a','MSC6b','MSC7a','MSC7b','MSC8a','MSC8b','MSC9a','MSC9b','MS10a','MSC10b'},'upper','omitselfloops');
% LWidths = abs((zscore(H.Edges.Weight/mean(H.Edges.Weight))));
% plot(H,'LineWidth',LWidths);
% title('Zscored edge weights')templ_match_dscalarwithassignments
if overlap ==0 % run overlap
    
    
    if exist('templ_match_dscalarwithassignments','var') == 1
        X{1} = diff_MIn_values;X{2} = MIn_templ;
        %X{1} = diff_MIn_values_info;X{2} = MIn_info;
        figure()
        options.shown_as='stairs';
        options.n_bins = 20;
        
        ct{1}='box';
        ct{2}='curve';
        ct{3}='stairs';
        ct{4}='contour';
        %options.n_bins=[];
        options.LineWidth=1.5;
        % options.n_bins=211;
        my_color=[27,158,119
            217,95,2
            117,112,179
            231,41,138]/255;
        clf
        
        %uncomment to plot both mono and dizygotic twins.
        %all_nulld = XD_muI{1};
        %all_nullm = XM_muI{1};
        %both_null = [all_nulld; all_nullm];
        %bothX{1} = both_null; bothX{2} = XM_muI{2}; bothX{3} = XD_muI{2};
        
        for i=1:4
            
            options.shown_as=ct{i};
            
            subplot(4,2,2*i-1)
            custom_hist(X,options)
            xlim([0.2 0.8])
            title (['shown as ' ct{i}])
            
            % providing your own color
            subplot(4,2,2*i)
            custom_hist(X,options,my_color)
            xlim([0.2 0.8])
            title (['shown as ' ct{i}])
        end
        %
        %
    else
    end
    
    
    if exist('infomap_dscalarwithassignments','var') == 1
        %X{1} = diff_MIn_values;X{2} = MIn_templ;
        X{1} = diff_MIn_values_info;X{2} = MIn_info;
        figure()
        options.shown_as='stairs';
        options.n_bins = 20;
        
        ct{1}='box';
        ct{2}='curve';
        ct{3}='stairs';
        ct{4}='contour';
        %options.n_bins=[];
        options.LineWidth=1.5;
        % options.n_bins=211;
        my_color=[27,158,119
            217,95,2
            117,112,179
            231,41,138]/255;
        clf
        
        %uncomment to plot both mono and dizygotic twins.
        %all_nulld = XD_muI{1};
        %all_nullm = XM_muI{1};
        %both_null = [all_nulld; all_nullm];
        %bothX{1} = both_null; bothX{2} = XM_muI{2}; bothX{3} = XD_muI{2};
        for i=1:4
            
            options.shown_as=ct{i};
            
            subplot(4,2,2*i-1)
            custom_hist(X,options)
            xlim([0.2 0.8])
            title (['shown as ' ct{i}])
            
            % providing your own color
            subplot(4,2,2*i)
            custom_hist(X,options,my_color)
            xlim([0.2 0.8])
            title (['shown as ' ct{i}])
        end
    end
    
else
    if exist('templ_match_dscalarwithassignments','var') == 1
        
        % all network colors
        all_colors = [
            
        255	0	0;
        0	0	153
        255	255	0
        0	255	0
        13	133	160
        50	50	50
        102	0	204
        102	255	255
        255	128	0
        178	102	255
        0	102	153
        102	255	102
        60	60	251
        200	200	200]/255;
    
    figure()
    for j=1: size(allscalars_templ,2)
        X{1} = diff_MIn_values(:,j);X{2} = MIn_templ(:,j);
        %X{1} = diff_MIn_values_info;X{2} = MIn_info;
        
        options.shown_as='stairs';
        options.n_bins = 30;
        
        ct{1}='box';
        ct{2}='curve';
        ct{3}='stairs';
        ct{4}='contour';
        %options.n_bins=[];
        options.LineWidth=1.5;
        % options.n_bins=211;
        %my_color=[27,158,119
        %    217,95,2
        %    117,112,179
        %    231,41,138]/255;
        null_color = [0 0 0];
        current_color = all_colors(j,:);
        my_color= [null_color; current_color];
        %my_color=all_colors(j);
        
        %uncomment to plot both mono and dizygotic twins.
        %all_nulld = XD_muI{1};
        %all_nullm = XM_muI{1};
        %both_null = [all_nulld; all_nullm];
        %bothX{1} = both_null; bothX{2} = XM_muI{2}; bothX{3} = XD_muI{2};
        
        options.shown_as=ct{1};
        
        % providing your own color
        subplot(4,4,j)
        custom_hist(X,options,my_color)
        xlim([0 0.9])
        title ([network_names{j}])
    end
    
    end
    
    
end
disp('Done running pairwise_mutalinfo.')
% %try
% % H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
% %
% % beta = 0 is a ring lattice, and beta = 1 is a random graph.
%
% % Connect each node to its K next and previous neighbors. This constructs
% % indices for a ring lattice.
% K = 10;
%
% s = repelem((1:Number_subjects)',1,K);
% t = s + repmat(1:K,Number_subjects,1);
% t = mod(t-1,Number_subjects)+1;
%
% % Rewire the target node of each edge with probability beta
% for source=1:Number_subjects
%     switchEdge = rand(K, 1) < beta;
%
%     newTargets = rand(Number_subjects, 1);
%     newTargets(source) = 0;
%     newTargets(s(t==source)) = 0;
%     newTargets(t(source, ~switchEdge)) = 0;
%
%     [~, ind] = sort(newTargets, 'descend');
%     t(source, switchEdge) = ind(1:nnz(switchEdge));
% end