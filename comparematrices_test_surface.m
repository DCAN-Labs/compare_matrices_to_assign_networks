function [new_subject_labels, output_cifti_scalar_name] = comparematrices_test_surface(input_cifti,output_cifti_name,method,data_type,cifti_enhancement,allow_overlap,overlap_method)
%10/25/18 Hermosillo R.

%This code uses a connectivity matrices, a template connectivity matrix, and label file, to try to assingn
%specific networks to the nodes of an individual.

%Currently this code only supports dconns pconns are in beta.

%what do you want to test?
close all
%data_type = 'dense'; % type dense or parcellated
%useeta =0; % calculates and eta value of every ROI to every ROI in the template. choose either correlation or eta
%usecorrelation =1; %runs a correlation of every ROI to every ROI in the template. choose either correlation or eta
%method = 'template_matching'; %select 'correlate' 'useeta' or 'template_matching'.
eta_method = 'average'; % eta_value
make_cifti_from_results =1; %output results as a viewable cifti.
%output_cifti_name = 'MSC01half1_to_template';
plot_data = 0; ROI = 50;
include_self = 0; % remove self from correlation (e.g. remove 1s from diagnoal when calculating correlation.
distribute_job =0; % Not supported. For use with dconns.  If distribute jobs =1, the code will attempt to split the 8.3bn comparisons (correlations or eta-squared) across 20 nodes. NOTE: Currently only tested with correlations.
%transform_data = 'Convert_FisherZ_to_r'; %cases:'Convert_FisherZ_to_r' or 'Convert_r_to_Pearons' or 'Convert_to_Zscores' %Used for template matching.  If your template values are in Pearson's r, but your data is Fisher Z transformed, this will convert your data from Fisher Z to Pearson. 
transform_data = 'Convert_to_Zscores';

%% Begin code.
if cifti_enhancement ==1
output_cifti_name = [output_cifti_name '_method_' method '_NE'];    
else
output_cifti_name = [output_cifti_name '_method_' method];
end
%input_cifti='/mnt/max/shared/projects/ADHD_comm_det/Infomap/concat_visits/8036-1_11_14_16_merged_SMOOTHED_2.55.dtseries.nii_all_frames_at_FD_0.2.dconn.nii';
%input_cifti='/mnt/max/shared/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/11480-1/20130813-SIEMENS_TrioTim-Nagel_K_Study/HCP_release_20161027/11480-1/MNINonLinear/Results/11480-1_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii';
%path_to_template_dconn = '/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii';


%addpath('/mnt/max/shared/code/external/utilities/Matlab_effect_size_toolbox/')
%wb_command = 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command'; % workbench command path
%addpath('/mnt/max/shared/code/external/utilities/Matlab_CIFTI')
%addpath(genpath('/mnt/max/shared/code/external/utilities/gifti-1.6'))
%addpath('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_scripts/')


%% Adding paths for this function
this_code = which('comparematrices_test_surface');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,10);

warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
wb_command=settings.path_wb_c; %path to wb_command
warning('on')

%% Start running code
% if usetestvector == 1
%     V1 = [ 1 1 2 1 2 3 3 1 1 2 7 7 7 7 7 7 9 9 9 9 9 9 9 ];
%     V2 = [ 2 2 4 2 4 1 4 2 2 3 7 7 7 7 7 7 9 9 9 9 9 9 9 ];
%     ncat = min(max(V1), max(template));
%     maxcategories = min(max(V1), max(template));
%
%     %accumarray(V1',1);
%     [~, lut1] = sort(accumarray(V1.', ones(size(V1))), 'descend') ;
%     mat_1networks = [1 1 1 1 2 2 3 3 4 4 4 4];
% end

if exist(input_cifti) == 0
    NOTE = ['cifti file does not exist']
    return
else
end

if strcmp(method,'template_matching')==1
else
    switch data_type
        
        case 'parcellated'
            disp('opening template cifti...')
            %cii=ciftiopen('/mnt/max/shared/projects/hcp_community_detection/Evan_test/Cifti_Community_Detection/distmat_creation/EUGEODistancematrix_XYZ_255interhem_unit8.pconn.nii',path_wb_c);
            %template_cii=ciftiopen('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii', wb_command);
            template_cii=ciftiopen(settings.path{4}, wb_command); %pconn path
            template_conn = template_cii.cdata;
            %template_cii=ciftiopen('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii', wb_command);
            template_cii=ciftiopen(settings.path{5}, wb_command); %parcellated networks labels path (pscalar).
            template_labels = template_cii.cdata;
            %additional template needed with more ROIs for saving (i.e. 352 vs 333).
            %saving_template =ciftiopen('/mnt/max/shared/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/10050-1/20100430-SIEMENS_TrioTim-Nagel_K_Study/HCP_release_20161027/10050-1/MNINonLinear/Results/10050-1_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii_to_Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii.pscalar.nii',wb_command);
            saving_template =ciftiopen(settings.path{7}, wb_command); % path to pscalar
            %subject_cii=ciftiopen('/mnt/max/shared/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/10050-2/20120418-SIEMENS-Nagel_K-Study/HCP_release_20161027/10050-2/MNINonLinear/Results/10050-2_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii', wb_command);
            subject_cii=ciftiopen(input_cifti, wb_command); %finally open the subject's data
            subject_conn = subject_cii.cdata;
            
            
        case 'dense'
            disp('opening template dconn...(1/4)')
            %cii=ciftiopen('/mnt/max/shared/projects/hcp_community_detection/Evan_test/Cifti_Community_Detection/distmat_creation/EUGEODistancematrix_XYZ_255interhem_unit8.pconn.nii',path_wb_c);
            %template_cii=ciftiopen('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii', wb_command);
            template_cii=ciftiopen(settings.path{16}, wb_command); %dconn path
            template_conn = single(template_cii.cdata);
            clear template_cii %save memory
            disp('opening dscalar of network labels...(2/4)');
            %template_cii=ciftiopen('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.dscalar.nii', wb_command);
            template_cii=ciftiopen(settings.path{16}, wb_command); %dense network labels (dscalar).
            template_labels = single(template_cii.cdata);
            num_surface_grey = length(template_labels);
            %additional template needed with more ROIs (i.e. 352 vs 333).
            disp('opening dscalar for saving...(3/4)');
            %saving_template =ciftiopen('/mnt/max/shared/code/internal/pipelines/HCP_release_20170910_v1.4/global/templates/91282_Greyordinates/91282_Greyordinates.dscalar.nii',wb_command);
            saving_template =ciftiopen(settings.path{15}, wb_command);
            %subject_cii=ciftiopen('/mnt/max/shared/data/study/ADHD/HCP/processed/ADHD_NoFMNoT2/10050-2/20120418-SIEMENS-Nagel_K-Study/HCP_release_20161027/10050-2/MNINonLinear/Results/10050-2_FNL_preproc_Gordon_subcortical.ptseries.nii_5_minutes_of_data_at_FD_0.2.pconn.nii', wb_command);
            disp('opening subject dconn...(4/4)');
            subject_cii=ciftiopen(input_cifti, wb_command); % finally opening subject's data
            subject_conn = single(subject_cii.cdata);
            clear subject_cii %save memory
            
        otherwise
            disp('Input type not support. Please select either parcellated or dense.')
    end
end

switch method
    case 'correlate'
        %% Correlation calculation (with self)
        correl = single(zeros(length(subject_conn),length(subject_conn))); %preallocate for speed
        if include_self ==1
            disp('calculating subject to template correlation for every ROI (~1 minute)');
            for subj_par_num=1:length(subject_conn)
                for templ_par_num=1:length(template_conn)
                    correl(subj_par_num,templ_par_num) = corr(subject_conn(:,subj_par_num),template_conn(:,templ_par_num)); %include self in correlation.
                    %correl(subject_parcel,template_parcel) = corr(subject_pconn([1:subject_parcel-1 subject_parcel+1:end],subject_parcel),template_pconn([1:template_parcel-1 template_parcel+1:end],template_parcel)); % do not include self in correlation
                    
                end
            end
            
            for i = 1:length(correl)
                %maximum_correl = max(correl(i,:));
                maximum_correl = max(correl(i,1:333)); %HARDCODE WARNING.  The code only looks for maxmium correlations among the parcellations for which there is a KNOWN group assingment.  Therefore it does not looks up the maximum correlations within subcortical regions.
                x(i) = find(correl(i,:) == maximum_correl);
            end
            
            new_subject_labels = template_labels(x); %get the labels based on the maximum corrlation values found above.
        else
            
            %% Correlation calulation (no self)
            disp('calculating subject to template correlation for every ROI (~1 minute for pconns.)');
            tic
            %correl = zeros(length(subject_conn),length(subject_conn));
            if distribute_job ==1
                cluster_env = parcluster();
                processingpool = parpool(cluster_env,numpools);
            end
            
            if distribute_job ==1
                %                 parfor subj_par_num=1:length(subject_conn)
                %                     for templ_par_num=1:length(template_conn)
                %
                %                         %correl(subject_parcel,template_parcel) = corr(subject_pconn(:,subject_parcel),template_pconn(:,template_parcel)); %include self in correlation.
                %                         if subj_par_num == templ_par_num
                %                             correl(subj_par_num,templ_par_num) = corr(subject_conn([1:subj_par_num-1 subj_par_num+1:end],subj_par_num),template_conn([1:templ_par_num-1 templ_par_num+1:end],templ_par_num)); % do not include self in correlation
                %                         elseif subj_par_num < templ_par_num
                %                             %removing the large fisher Z value (e.g. 7.2) from each vector in the correlation results in an offset in the vectors when the ROI number is not equal.
                %                             % To handle the removal of outliers from each vector, while
                %                             % maintaining the paired association of the ROIs, we remove
                %                             % both the outlier value and its corredponding ROI in the
                %                             % other vector. (i.e. instead of the correlation performed
                %                             % 352 "X" values and 352 "Y" values, we remove the outliers
                %                             % ROI rfom both so the correlation is between 350 X values nad
                %                             % 350 Y values.
                %                             correl(subj_par_num,templ_par_num) = corr(subject_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],subj_par_num),template_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],templ_par_num)); % do not include self or other in correlation
                %                         elseif subj_par_num > templ_par_num
                %                             correl(subj_par_num,templ_par_num) = corr(subject_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],subj_par_num),template_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],templ_par_num)); % do not include self or other in correlation
                %
                %                         else
                %                             disp('something went wrong with the number of ROI calculation')
                %                         end
                %
                %                     end
                %                 end
                disp(' Distributing job not currently supported');
                return
                toc
                
                %close the pool all other operations only need one core after
                if parallel_processing
                    delete(processingpool);
                end
                
            else %do not distribute job
                for subj_par_num=1:length(subject_conn)
                    for templ_par_num=1:length(template_conn)
                        
                        %correl(subject_parcel,template_parcel) = corr(subject_pconn(:,subject_parcel),template_pconn(:,template_parcel)); %include self in correlation.
                        if subj_par_num == templ_par_num
                            correl(subj_par_num,templ_par_num) = corr(subject_conn([1:subj_par_num-1 subj_par_num+1:end],subj_par_num),template_conn([1:templ_par_num-1 templ_par_num+1:end],templ_par_num)); % do not include self in correlation
                        elseif subj_par_num < templ_par_num
                            %removing the large fisher Z value (e.g. 7.2) from each vector in the correlation results in an offset in the vectors when the ROI number is not equal.
                            % To handle the removal of outliers from each vector, while
                            % maintaining the paired association of the ROIs, we remove
                            % both the outlier value and its corredponding ROI in the
                            % other vector. (i.e. instead of the correlation performed
                            % 352 "X" values and 352 "Y" values, we remove the outliers
                            % ROI rfom both so the correlation is between 350 X values nad
                            % 350 Y values.
                            correl(subj_par_num,templ_par_num) = corr(subject_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],subj_par_num),template_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],templ_par_num)); % do not include self or other in correlation
                        elseif subj_par_num > templ_par_num
                            correl(subj_par_num,templ_par_num) = corr(subject_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],subj_par_num),template_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],templ_par_num)); % do not include self or other in correlation
                            
                        else
                            disp('something went wrong with the number of ROI calculation')
                        end
                        if templ_par_num ==1
                            if subj_par_num==1 || subj_par_num==10000 || subj_par_num==20000 || subj_par_num==30000 || subj_par_num==40000 || subj_par_num==50000 || subj_par_num==60000 || subj_par_num==70000 || subj_par_num==80000 || subj_par_num==90000
                                disp(subj_par_num);toc;
                            end
                        end
                        
                    end
                end
                toc
            end
            
            %% Find maximum correlation in template
            disp('finding maximum  correlation');
            switch data_type
                case 'parcellated'
                    for i = 1:length(correl)
                        %maximum_correl = max(correl(i,:));
                        maximum_correl = max(correl(i,1:333)); %HARDCODE WARNING.  The code only looks for maxmium correlations among the parcellations for which there is a KNOWN group assingment.  Therefore it does not looks up the maximum correlations within subcortical regions.
                        x(i) = find(correl(i,:) == maximum_correl);
                    end
                case 'dense'
                    for i = 1:length(correl)
                        %maximum_correl = max(correl(i,:));
                        maximum_correl = max(correl(i,1:num_surface_grey)); %HARDCODE WARNING.  The code only looks for maxmium correlations among the parcellations for which there is a KNOWN group assingment.  Therefore it does not looks up the maximum correlations within subcortical regions.
                        x(i) = find(correl(i,:) == maximum_correl);
                    end
                otherwise
                    disp('Error: data type either not supported or not found.  Please select parcellated or dense.');
            end
            
            new_subject_labels = template_labels(x); %get the labels based on the maximum corrlation values found above.
            
            %figure()
            if plot_data ==1
                figure()
                subplot(2,5,1)
                imagesc(template_conn)
                xlabel('template ROI'); ylabel('template ROI');
                subplot(2,5,6)
                imagesc(subject_conn)
                xlabel('subject ROI'); ylabel('subject ROI');
                imagesc(correl)
                ylabel('subject ROI including self'); xlabel('template ROI including self');
                subplot(2,5,3)
                scatter(template_labels(1:333),new_subject_labels(1:333));
                ylabel('subject ROI labels self'); xlabel( 'template ROI labels including self');
                subplot(2,5,4)
                scatter(subject_conn(:,ROI_num),template_conn(:,ROI_num))
                ylabel('subject ROI(example)'); xlabel('template ROI example');
                subplot(2,5,5)
                histogram(template_labels,'FaceAlpha',0.5)
                hold on
                histogram(new_subject_labels,'FaceAlpha',0.5)
                ylabel('number of parcellations'); xlabel('network number');
                subplot(2,5,7)
                imagesc(correl)
                ylabel('subject ROI excluding self'); xlabel('template ROI excluding self');
                subplot(2,5,8)
                scatter(template_labels(1:333),new_subject_labels(1:333));
                ylabel('subject ROI labels excluding self'); xlabel( 'template ROI labels exclusing self');
                subplot(2,5,9)
                scatter(subject_conn([1:ROI_num-1 ROI_num+1:end],ROI),template_conn([1:ROI_num-1 ROI_num+1:end],ROI_num))
                ylabel('subject ROI(example)'); xlabel('template ROI example');
                subplot(2,5,10)
                histogram(template_labels,'FaceAlpha',0.6)
                hold on
                histogram(new_subject_labels,'FaceAlpha',0.6)
                ylabel('number of parcellations'); xlabel('network number');
            end
        end
        
        
        %% Eta squared calulation (no self)
    case 'useeta'
        %if useeta ==1
        tic
        try
            disp('Eta matrix file found. Loading.')
            load('/mnt/max/shared/code/external/utilities/Matlab_effect_size_toolbox/11480-1_to_template.mat');
        catch
            disp('Eta matrix file not found.')
        end
        
        if exist('eta_mat') == 0
            disp('calculating subject to template eta squared for every ROI (~1hr)');
            tic
            eta_mat = zeros(length(subject_conn),length(subject_conn));
            for subj_par_num=1:length(subject_conn)
                for templ_par_num=1:length(template_conn)
                    %correl(subject_parcel,template_parcel) = corr(subject_pconn(:,subject_parcel),template_pconn(:,template_parcel)); %include self in correlation.
                    if subj_par_num == templ_par_num
                        paired_eta_mat = [subject_conn([1:subj_par_num-1 subj_par_num+1:end],subj_par_num) template_conn([1:templ_par_num-1 templ_par_num+1:end],templ_par_num)]; % do not include self in correlation
                    elseif subj_par_num < templ_par_num
                        %removing the large fisher Z value (e.g. 7.2) from each vector in the correlation results in an offset in the vectors when the ROI number is not equal.
                        % To handle the removal of outliers from each vector, while
                        % maintaining the paired association of the ROIs, we remove
                        % both the outlier value and its corredponding ROI in the
                        % other vector. (i.e. instead of the correlation performed
                        % 352 "X" values and 352 "Y" values, we remove the outliers
                        % ROI rfom both so the correlation is between 350 X values nad
                        % 350 Y values.
                        paired_eta_mat = [subject_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],subj_par_num) template_conn([1:subj_par_num-1 subj_par_num+1:templ_par_num-1  templ_par_num+1:end],templ_par_num)]; % do not include self or other in correlation
                    elseif subj_par_num > templ_par_num
                        paired_eta_mat = [subject_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],subj_par_num) template_conn([1:templ_par_num-1 templ_par_num+1:subj_par_num-1  subj_par_num+1:end],templ_par_num)]; % do not include self or other in correlation
                        
                    else
                        disp('something went wrong with the number of ROI calculation')
                    end
                    [paired_results,~] = mes1way(paired_eta_mat,'eta2');
                    eta_mat(subj_par_num,templ_par_num)=paired_results.eta2;
                    if templ_par_num ==1
                        if subj_par_num==1 || subj_par_num==10000 || subj_par_num==20000 || subj_par_num==30000 || subj_par_num==40000 || subj_par_num==50000 || subj_par_num==60000 || subj_par_num==70000 || subj_par_num==80000 || subj_par_num==90000
                            disp(subj_par_num);toc;
                        end
                    end
                end
                disp(subj_par_num);
            end
            toc
        else%exist(eta_mat)
        end %exist(eta_mat)
        
        
        %group eta values by network.
        num_networks = max(template_labels);
        for i =1:num_networks
            network_indices{i} = find(template_labels == i);
        end
        eta_network_summary = zeros(length(eta_mat),num_networks);
        
        switch eta_method
            
            case 'average'
                disp('Using "average" as specified by user.')
                for j = 1:length(eta_mat)
                    for k = 1:num_networks
                        network_eta = eta_mat(j,network_indices{1,k});
                        eta_network_summary(j,k) = mean(network_eta);
                    end
                end
                
            case 'median'
                disp('Using "median" as specified by user.')
                for j = 1:length(eta_mat)
                    for k = 1:num_networks
                        network_eta = eta_mat(j,network_indices{1,k});
                        eta_network_summary(j,k) = mean(network_eta);
                    end
                end
                
            case 'max'
                disp('Using "max" as specified by user.')
                for j = 1:length(eta_mat)
                    for k = 1:num_networks
                        network_eta = eta_mat(j,network_indices{1,k});
                        eta_network_summary(j,k) = max(network_eta);
                    end
                end
            otherwise
                error('Error: method is no found or allowed')
        end
        
        for i = 1:length(eta_mat)
            %maximum_correl = max(correl(i,:));
            %network_eta = max(eta_mat(i,1:333)); %HARDCODE WARNING.  The code only looks for maxmium correlations among the parcellations for which there is a KNOWN group assingment.  Therefore it does not looks up the maximum correlations within subcortical regions.
            maximum_eta = max(eta_network_summary(i,1:num_networks)); %HARDCODE WARNING.  The code only looks for maxmium correlations among the parcellations for which there is a KNOWN group assingment.  Therefore it does not looks up the maximum correlations within subcortical regions.
            new_subject_labels(i) = find(eta_network_summary(i,1:num_networks) == maximum_eta);
            new_subject_labels = new_subject_labels';
        end
        %else%useeta
        %end %useeta
        
    case 'template_matching'
        disp('Template matching method selected. Running...');
        if exist([output_cifti_name '.mat']) == 0 % if the matrix already exists skip
            disp([output_cifti_name '.mat file (contains Eta_to_template_vox) Eta_to_template_vox .mat file not found for this subject found. Running template matching (~1-2 hrs).'])
            %networks_template_path = settings.path{12};
            networks_template_path = settings.path_template_nets;
            networks_template_path = ('/home/exacloud/lustre1/fnl_lab/projects/INFANT/GEN_INFANT/Luci_template/infant_template/seedmaps_UCI_smoothed_dtseries_all_networksZscored.mat');
            cifti_output_folder = pwd;
            [new_subject_labels,output_cifti_scalar_name] = template_matching_LM(input_cifti,data_type, networks_template_path,transform_data,output_cifti_name,cifti_output_folder,wb_command,make_cifti_from_results,allow_overlap,overlap_method);
            %[eta_to_template_vox, eta_subject_index,output_cifti_scalar_name] = template_matching_RH(subjectlist, data_type, template_path,transform_data,output_cifti_name, cifti_output_folder, wb_command, make_cifti_from_results)

            %new_subject_labels = eta_subject_index;
            %save([output_cifti_name '.mat'],'eta_to_template_vox','eta_subject_index','-v7.3')
            %saving_template =ciftiopen(settings.path{16}, wb_command);   
        else
            disp([output_cifti_name '.mat file (contains Eta_to_template_vox) found for this subject found. Template matching has alreday been run for subject. Loading data.'])
            %networks_template_path = settings.path{12};
            networks_template_path = settings.path_template_nets;
            saving_template =ciftiopen(settings.path{16}, wb_command);
            load([output_cifti_name '.mat'])  
            new_subject_labels = eta_subject_index;
            if make_cifti_from_results == 1
                disp('saving file to cifti')
                saving_template.cdata = single(new_subject_labels);
                %addpath('/mnt/max/shared/code/internal/utilities/corr_pt_dt/support_files');
                disp('Saving new scalar')
                save(saving_template, [output_cifti_name '.gii'], 'ExternalFileBinary')
                disp('Converting Zscored scalar .gii to .nii')
                switch data_type
                    case 'parcellated'
                        output_cifti_scalar_name  = [output_cifti_name '.pscalar.nii' ];
                    case 'dense'
                        output_cifti_scalar_name  = [output_cifti_name '.dscalar.nii' ];
                end
                unix([wb_command ' -cifti-convert -from-gifti-ext ' output_cifti_name '.gii ' output_cifti_scalar_name ]);
                disp('Removing .gii')
                unix(['rm -f ' output_cifti_name '.gii']);
                unix(['rm -f ' output_cifti_name '.dat']);
            else
            end
        end
    otherwise
        Disp('"method_type" either not supported or not found.  Please select "average" "useeta" or "template_matching"  (in single quotes)');
end



disp('done')
