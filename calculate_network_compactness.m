function [scores,cluster_num_pvalue, compactness_pvalue,avg_compact_pval,avg_num_clusters_pval] = calculate_network_compactness(network_suface_area_mat_file,network_perimeter_mat_file,method,run_get_pairwise_stats)
% This function calcualtes the compactness of a network if given a list of surface areas and networks for each subject.

%inputs are:
% 'network_suface_area_mat_file' = a .mat file of surface areas in mm^2
% 'network_perimeter_mat_file' = a .mat file of network perimeter (i.e. border length) in mm.
% 'method' = and desired method type.  See Below.
% 'run_get_pairwise_stats' = get pairwise p-values for compactness and number of clusters.

%outputs are:
%score = compactness score
%if run_get_pairwise_stats is set to true (==1), then the code will coade
%can calculate the fllowing statistics:
%cluster_num_pvalue = test if the number of clusters of each network is different between the even and odd subjects.
%compactness_pvalue = test if the compactness each network is different between the even and odd subjects.
%avg_compact_pval = test if the compactness across all network is different between the even and odd subjects.
%avg_num_clusters_pval = test if the average number of clusters across all network is different between the even and odd subjects.

% The term Gerrymandering refers to the act of manipulating the boundaries
% of voting districts to achieve some political advantage. The term was
% coined during tenure Massachusetts Governor Elbridge Gerry, who in 1812
% redrew the voting districts for the Massachusetts State Senate to favor
% his own party. One district caught the attention of the Boston Gazette,
% who published a political cartoon likening the district�s shape to that
% of a salamander and labeling the phenomenon �The Gerry-mander� after the
% Governor.

%six of the most frequently used measures of compactness used by academic
% researchers: (1) Polsby-Popper (Polsby and Popper, 1991); (2) Schwartzberg
% (1965); (3) Reock (1961); (4) Convex Hull; (5) X-Symmetry; and (6)
% Length-Width Ratio (C.C. Harris, 1964).

%Polsby-Popper
%The Polsby-Popper (PP) measure (polsby & Popper, 1991) is the ratio of the area of the district (AD) to the area of a circle whose circumference is equal to the perimeter of the district (PD). A district�s Polsby-Popper score falls with the range of [0,1] and a score closer to 1 indicates a more compact district.
%PP = 4pi x (A/(P^2))

%Schwartzberg
%The Schwartzberg score (S) compactness score is the ratio of the perimeter of the district (PD) to the circumference of a circle whose area is equal to the area of the district. A district�s Schwartzberg score as calculated below falls with the range of [0,1] and a score closer to 1 indicates a more compact district.
%S = 1/P/C = 1/(P/(2pi*sqrt((A/pi))));

%Reock Score
%The Reock Score (R) is the ratio of the area of the district AD to the area of a minimum bounding cirle (AMBC) that encloses the district�s geometry. A district�s Reock score falls within the range of [0,1] and a score closer to 1 indicates a more compact district.
%R = ((A)/(AMBC)); NOTE: not currently supported as the calculation of the  minimum bounding circle
%requires that that points lie on a 2D plane.

%Start
load(network_suface_area_mat_file,'network_surfarea'); %
load(network_perimeter_mat_file,'network_lengths_for_each_sub');

scores = zeros(size(network_surfarea,1), size(network_surfarea,2)-2);

for i = 1:size(network_surfarea,1)
    disp(i)
    this_subs_surface_areas = network_surfarea(i,:);
    %this_subs_perimeters = network_lengths_for_each_sub(i,:);
    
    %remove non-existing networks #4 and 6
    this_subs_surface_areas(:,4)  = [];
    this_subs_surface_areas(:,5)  = [];
    
    for j = 1:size(this_subs_surface_areas,2)
                        total_length_this_sub = network_lengths_for_each_sub{i, 1}{j, 4}+ network_lengths_for_each_sub{i, 2}{j, 4};

        switch method
            
            case 'polsbypopper'
                %disp('method is polsbypopper')
                %combine hemisphere lengths
                scores(i,j) = 4*pi*(this_subs_surface_areas(1,j)/(total_length_this_sub^2));
                
            case 'schwartzberg'
                 scores(i,j) = 1/(total_length_this_sub/(2*pi*(sqrt(this_subs_surface_areas(1,j)/pi))));
               
            otherwise
                disp('method is not supported')
        end
    end
end
if run_get_pairwise_stats ==1
    [cluster_num_pvalue, compactness_pvalue,avg_compact_pval,avg_num_clusters_pval] = get_pairwise_stats(scores,network_lengths_for_each_sub);
else
    cluster_num_pvalue = []; 
    compactness_pvalue=[];
    avg_compact_pval=[];
    avg_num_clusters_pval=[];
end
end

function [cluster_num_pvalue, compactness_pvalue,avg_compact_pval,avg_num_clusters_pval]=get_pairwise_stats(scores,network_lengths_for_each_sub)
disp('Testing each element of the matrix, pairwise')
twin1_indices = 1:2:size(scores,1); % build indices vector
twin2_indices = 2:2:size(scores,1); % build indices vector

%preallocate
twin1_scores =zeros(size(scores,1),size(scores,2));
twin2_scores =zeros(size(scores,1),size(scores,2));
twin1_cluster_num =zeros(size(scores,1),size(scores,2));
twin2_cluster_num =zeros(size(scores,1),size(scores,2));
compactness_pvalue =zeros(size(twin1_indices,2),1);
cluster_num_pvalue =zeros(size(twin1_indices,2),1);

k=1;
for  i =  twin1_indices
    twin1_scores(k,:) =  scores(i,:);
    twin1_cluster_num(k,:) = (cell2mat(network_lengths_for_each_sub{i, 1}(:,2)) + cell2mat(network_lengths_for_each_sub{i, 2}(:,2)))'; % combine number of clusters from each hemisphere.
    k=k+1;
end
k=1;
for  i =  twin2_indices
    twin2_scores(k,:) =  scores(i,:);
    twin2_cluster_num(k,:) = (cell2mat(network_lengths_for_each_sub{i, 1}(:,2)) + cell2mat(network_lengths_for_each_sub{i, 2}(:,2)))'; % combine number of clusters from each hemisphere.
    k=k+1;
end

for j = 1:size(twin1_indices,2) % go through each network
    [~, compact_p, ~, ~]=ttest(twin1_scores(:,j),twin2_scores(:,j));
    [~, cluster_num_p, ~, ~]=ttest(twin1_cluster_num(:,j),twin2_cluster_num(:,j));
    %[hypoth, pval, CI, tstats]=ttest(twin1_scores,twin2_scores);
    %paired_hypoth(j,k)=hypoth;
    compactness_pvalue(j,1)=compact_p;
    cluster_num_pvalue(j,1)=cluster_num_p;
    %pairedCI{j,k}= CI;
    %paired_t_test{j,k}= tstats;
end

avg_compact1=mean(twin1_scores,2); %average network compactness across networks
avg_num1=mean(twin1_cluster_num,2);
avg_compact2=mean(twin2_scores,2);
avg_num2=mean(twin2_cluster_num,2);

[~, avg_compact_pval, ~, ~]=ttest(avg_compact1,avg_compact2);
[~,avg_num_clusters_pval, ~, ~]=ttest(avg_num1,avg_num2);

end