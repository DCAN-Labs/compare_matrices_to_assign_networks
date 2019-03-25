%load('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/surface_area/lucianov_test_cleaned.mat');%load surface and volume data
%load('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/surface_area/ADHDsymptoms.mat')
%dscalarwithassignments = importdata('/mnt/max/shared/projects/NIGGTWINS/WTO/Experiments/Template_matching/template_matching_dscalars/template_matching_cleaned_dscalar.conc');
dscalarwithassignments = importdata('/mnt/max/shared/projects/midnight_scan_club/template_matching/bothhalvesdsclars.conc');

%% Add necessary paths
addpath ('/mnt/max/shared/code/internal/analyses/compare_matrices')
this_code = which('surfaceareafromgreyordinates');
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

%Lmidthicknessfile,Rmidthicknessfile,dscalarwithassignments

%conc = strsplit(dscalarwithassignments, '.');
%conc = char(conc(end));
%if strcmp('conc',conc) == 1
%    dscalarwithassignments = importdata(dscalarwithassignments);
%else
%    dscalarwithassignments = {dscalarwithassignments};
%end

%check to make sure that surface files exist
for i = 1:length(dscalarwithassignments)
    if exist(dscalarwithassignments{i}) == 0
        NOTE = ['Subject dscalar ' num2str(i) ' does not exist']
        disp(dscalarwithassignments{i});
        return
    else
    end
end
disp('All dscalars exist continuing ...')

Number_subjects = length(dscalarwithassignments);
Number_network_surfarea = 16;

for i = 1:Number_subjects
allscalars_temp = ciftiopen(dscalarwithassignments{i},wb_command);
allscalars(:,i) = single(allscalars_temp.cdata);
end

% for i = 1:2:26 
%     subject1 = network_surfarea (i,:);
%     subject2 = network_surfarea (i+1,:);
%     network_difference_surfarea(round(i/2),:) = subject1 - subject2;
% end
% 
% 
% Number_subjects = 26;
% Number_network_surfarea = 16;
% for i = 1:2:26 
%     subject1 = network_volume (i,:);
%     subject2 = network_volume (i+1,:);
%     network_difference_volume(round(i/2),:) = subject1 - subject2;
% end

% all_network_difference_surfarea = zeros(26,26,16);
% 
 twin1_indices = 1:2:length(dscalarwithassignments);
 twin2_indices = 2:2:length(dscalarwithassignments);

% for j = 1:16 %go through every networks
%     for i = 1:26 %go through every subject
%         for k = 1:26 %go through every other subject
%             if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
%                 %if i < k
%                 subject1 = network_surfarea (i,j); eta_net_assign1 = allscalars(i);
%                 subjectk = network_surfarea (k,j); eta_net_assignk = allscalars(k);
%                 all_network_difference_surfarea(i,k,j) = abs(subject1 - subjectk);
%             end
%         end
%     end
% end


% for i=1:16   %go through every network
% diffmatrix = all_network_difference_surfarea(:,:,i);
% diff_indices = find(diffmatrix~=0);
% diff_values = diffmatrix(diff_indices);
% a = mean(network_difference_surfarea(:,i));
% proportion_of_diff_below_paired_average_diff(i) = length(find(diffmatrix(diff_indices)<a))/length(diff_indices);
% %histogram(diff_values,30); hold on; histogram(network_difference_surfarea(:,i),30,'FaceColor','r');
% %hist(diffmatrix(diff_indices),30); hold on; bar(a,16,50,'r')
% end

%calculate mutual information for real pairs
for i = 1:2:Number_subjects %number of subjects
    %Step 3: calculate mutualinformation
    %idx = 1:2:26;
    %muI(i,1) = MutualInformation(eta_net_assign{i},eta_net_assign{i+1}); %Mutual information
    %[VIn(i,1), MIn(i,1)] = partition_distance(eta_net_assign{i},eta_net_assign{i+1}); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
    muI(round(i/2),1) = MutualInformation(allscalars(:,i),allscalars(:,i+1)); %Mutual information
    [VIn(round(i/2),1), MIn(round(i/2),1)] = partition_distance(allscalars(:,i),allscalars(:,i+1)); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)

end

%calculate mututal informatin to all possible combinations, except pairs
for i = 1:Number_subjects %go through every subject
    for k = 1:Number_subjects %go through every other subject
        if i < k && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
            eta_net_assign1 = allscalars(:,i);
            eta_net_assignk = allscalars(:,k);
                unpaired_muI(i,k) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                [unpaired_VIn(i,k), unpaired_MIn(i,k)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
        end
    end
end

%calculate mututal informatin to all possible combinations, including pairs
for i = 1:Number_subjects %go through every subject
    for k = 1:Number_subjects %go through every other subject
        %if i < k% && twin1_indices(round(i/2))+1 ~=  twin2_indices(round(k/2)) %only calculate differences for non-twins.
            eta_net_assign1 = allscalars(:,i);
            eta_net_assignk = allscalars(:,k);
                allposs_muI_pairs(i,k) = MutualInformation(eta_net_assign1,eta_net_assignk); %Mutual information
                [allposs_w_pairs_VIn(i,k), allposs_w_pairs_MIn(i,k)] = partition_distance(eta_net_assign1,eta_net_assignk); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix)
       %end
    end
end


diff_muI_matrix = unpaired_muI;
%[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
diff_muI_indices = find(unpaired_muI~=0);
diff_muI_values = diff_muI_matrix(diff_muI_indices);
%a = mean(abs(muI));b = mean(diff_muI_values); c = std(diff_muI_values);
%proportion_of_mu_I_diff_below_paired_average_diff = length(find(diff_muI_values<a))/length(diff_muI_indices);

diff_VIn_matrix = unpaired_VIn;
%[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
diff_VIn_indices = find(unpaired_VIn~=0);
diff_VIn_values = diff_VIn_matrix(diff_VIn_indices);
%a = mean(abs(VIn));b = mean(diff_VIn_values); c = std(diff_VIn_values);
%proportion_of_mu_I_diff_below_paired_average_diff = length(find(diff_VIn_values<a))/length(diff_VIn_indices);

diff_MIn_matrix = unpaired_MIn;
%[diff_muI_indices_row,diff_muI_indices_colum] = find(allposs_muI~=0);
diff_MIn_indices = find(unpaired_MIn~=0);
diff_MIn_values = diff_MIn_matrix(diff_MIn_indices);
%a = mean(abs(MIn));b = mean(diff_MIn_values); c = std(diff_MIn_values);
%proportion_of_mu_I_diff_below_paired_average_diff =
%length(find(diff_MIn_values<a))/length(diff_MIn_indices);


%t = (a - b) / (c/sqrt(length(diff_muI_indices)));
%hist(diff_muI_values,30); hold on; bar(a,16,'BarWidth',0.01,'FaceColor','r');
figure(); histogram(diff_muI_values,30); hold on; histogram(muI,30,'FaceColor','r');xlabel('Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast')
figure(); histogram(diff_VIn_values,30); hold on; histogram(VIn,30,'FaceColor','r');xlabel('Variation of Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast')
figure(); histogram(diff_MIn_values,30); hold on; histogram(MIn,30,'FaceColor','r');xlabel('Normalized Mutual Information (bits)');ylabel('count');legend({'All possible half pairs','Real half pairs'},'Location','northeast')
G = graph(allposs_muI_pairs,{'MSC1a',',MSC1b','MSC2a',',MSC2b','MSC3a',',MSC3b','MSC4a',',MSC4b','MSC5a',',MSC5b','MSC6a',',MSC6b','MSC7a',',MSC7b','MSC8a',',MSC8b','MSC9a',',MSC9b','MS10a',',MSC10b'},'upper','omitselfloops');
figure(); plot(G);
disp('Done running pairwise_mutalinfo.')

figure()
plot(muI)
plot(muI,'MarkerSize',10,'MarkerEdgeColor','b')
boxplot(info_muI')


%try
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
K = 10;

s = repelem((1:Number_subjects)',1,K);
t = s + repmat(1:K,Number_subjects,1);
t = mod(t-1,Number_subjects)+1;

% Rewire the target node of each edge with probability beta
for source=1:Number_subjects    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(Number_subjects, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);


%catch
%disp(could not build Watts-Strogatz model')
%end


