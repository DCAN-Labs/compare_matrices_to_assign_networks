function [ proportion_time, transition_matrix,proportions_matrix] = network_state(network_weights_filepath, dtseries_file_path,motion_mask,output_mat_name,image_limit)

%R. Hermosillo 07/14/2023
% This function works by comparing each frame of the time series to a series of templates.
% These can be the template that are used for template matching, 
% or they can simply be a given subject's weights from template matching.


% %Inputs are
% 1) nework weights (string): the path to the .mat file containing the weights from template matching.
% 2) dtseries_file_path (string): The path to the time series to which you want to identfy the framewise pairing to the template.
% 3) motion_mask (string): the path to the motion .txt  file that contains the motion censor of the frames to use. (You can get this from cifti_conn_matrix).
% 4) output_mat_name(string): some putput name excluding the '.mat' part.
% 
% update 01/17/ Added output file name so that you can export the proportion_time and save the results.


%harcodes for testing
%load(network_weights_filepath,'eta_to_template_vox');
run_locally=0;
only_use_positive_data=1;
include_residuals=0;
do_spline_interpolation_on_nets=0;
do_spline_interpolation_on_dtseries=1;
paddingLength =5;
do_motion_censoring =1;
use_cortex_only =1;
k=15;

%some hardcodes
%load('/panfs/jay/groups/14/znahas/shared/projects/FEAST/precision_networks/merged_sessions/sub-FEAST0003_ses-merged_Zscored.mat','eta_to_template_vox');
%load('C:\Users\hermosir\Documents\test_ciftis\VA_METH_TAAR_dtseries\sub-TAAR711_ses-combined_merged_timeseries_template_matched_Zscored.mat','eta_to_template_vox');
load(network_weights_filepath);

%dtseries_file_path='/panfs/jay/groups/14/znahas/shared/projects/FEAST/derivatives/06292023/sub-FEAST0003_ses-01/xcp_d/sub-FEAST0003/ses-01/func/sub-FEAST0003_ses-01_task-restMENORDICrmnoisevols_space-fsLR_den-91k_desc-interpolated_bold_spatially_interpolated.dtseries.nii';
%dtseries_file_path='/panfs/jay/groups/14/znahas/shared/projects/FEAST/dconns/dconns_no_outlier_detection/sub-FEAST0003_ses-01_task-restMENORDICrmnoisevols_space-fsLR_den-91k_desc-interpolated_bold_spatially_interpolated_SMOOTHED_2.55.dtseries.nii';
%dtseries_file_path='/home/znahas/shared/projects/FEAST/dconns/separate_sessions/sub-FEAST0003_ses-01_task-restMENORDICrmnoisevols_space-fsLR_den-91k_desc-interpolated_bold_spatially_interpolated_SMOOTHED_2.55.dtseries.nii';
%dtseries_file_path='C:\Users\hermosir\Documents\test_ciftis\VA_METH_TAAR_dtseries\sub-TAAR711_ses-combined_task-restV3B_bold_desc-filtered_timeseries_spatially_interpolated.dtseries.nii';

%motion_mask=importdata('/panfs/jay/groups/14/znahas/shared/projects/FEAST/dconns/dconns_no_outlier_detection/sub-FEAST0003_ses-01_task-restMENORDICrmnoisevols_desc-dcan_qc_power_2014_FD_only.mat_0.2_cifti_censor_FD_vector_All_Good_Frames.txt');
%motion_mask=importdata('/home/znahas/shared/projects/FEAST/dconns/separate_sessions/sub-FEAST0003_ses-01_task-restMENORDICrmnoisevols_desc-dcan_qc_power_2014_FD_only.mat_0.2_cifti_censor_FD_vector_All_Good_Frames.txt');
%motion_mask=importdata('C:\Users\hermosir\Documents\test_ciftis\VA_METH_TAAR_dtseries/sub-TAAR711_ses-combined_task-restV3B_bold_desc-filtered_motion_mask_outliers.mat');
if strcmp(motion_mask(end-3:end),'.mat') ==1
    motion_mask=load(motion_mask);
    motion_mask=logical(motion_mask{1,21}.combined_removal);
    motion_mask=~motion_mask; %invert the motion mask %one equals remove
else
    motion_mask=importdata(motion_mask);
    motion_mask=logical(motion_mask); % ones equal keep.
end
network_names = {   'DMN'    'Vis'    'FP'  'nonet1'  'DAN'  'nonet2'    'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'  'nonet3' 'SCAN'};
eta_to_template_vox=double(eta_to_template_vox);

if use_cortex_only ==1
    eta_to_template_vox =  eta_to_template_vox(1:59412,1:end);
end
eta_to_template_vox(:,4)=[];network_names(4)=[];
eta_to_template_vox(:,5)=[];network_names(5)=[];
eta_to_template_vox(:,15)=[];network_names(15)=[];

%add cifti paths
if run_locally ==1
    %Some hardcodes:
    wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
    %support_folder='C:\Users\hermosir\Documents\repos\support_folder';
else
    this_code = which('template_matching_RH');
    [code_dir,~] = fileparts(this_code);
    support_folder=[code_dir '/support_files']; %find support files in the code directory.
    addpath(genpath(support_folder));
    settings=settings_comparematrices;%
    np=size(settings.path,2);
    disp('Attempting to add neccesaary paths and functions.')
    warning('off') %supress addpath warnings to nonfolders.
    for i=1:np
        addpath(genpath(settings.path{i}));
    end
    warning('on')
    % Check if wb_command has been provided
    if ~exist('wb_command', 'var') || isempty(wb_command)
        % If wb_command is not provided or is empty, set the default path
        wb_command = settings.path_wb_c; %path to wb_command
    end
    addpath(genpath('/projects/standard/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti/'));
end

%% Step 1 - Open timseries file
%[filepath,filenamewext1,ext2] = fileparts(dtseries_file_path);
%[~,filename,ext1] = fileparts(filenamewext1);
disp('Loading cifti...')
cii = ciftiopen(dtseries_file_path,wb_command);

timeseries = cii.cdata;
if use_cortex_only ==1
    timeseries=timeseries(1:59412,1:end);
end

num_time_frames = size(timeseries, 2);
num_networks = size(eta_to_template_vox,2);

if do_motion_censoring ==1
    motion_mask=logical(motion_mask);
else
    motion_mask=ones(num_time_frames,1);
    
end


if include_residuals ==1
    proportions_matrix = zeros(num_time_frames,num_networks +1); % add an extra colum for the residuals.
else
    proportions_matrix = zeros(num_time_frames,num_networks); % add an extra colum for the residuals.
end

%%  Do interpolation on dtseries
if do_spline_interpolation_on_dtseries==1
    % Define time vector
    t = 1:size(timeseries, 2)+2*(paddingLength);
    % Zero-pad the data to anchor the endpoints
    paddingLength = 5; % Example: Add 5 zeros on each side
    paddedTimeseries = [zeros(size(timeseries, 1), paddingLength), timeseries, zeros(size(timeseries, 1), paddingLength)];
    % Update motion mask to account for padding
    motion_mask_padded = [ones(paddingLength, 1); motion_mask; ones(paddingLength, 1)];
    validIndices = motion_mask_padded == 1;
    
    for i= 1:size(timeseries,2)+(2*paddingLength)
        if motion_mask_padded(i) == 0
            paddedTimeseries(:,i)=zeros(size(timeseries,1),1);
        else
        end
    end
    
    % Perform spline interpolation row by row
    interpolatedTimeseries = paddedTimeseries; % Copy padded data
    for row = 1:size(timeseries, 1)
        validTime = t(validIndices); % Time points corresponding to valid data
        validData = paddedTimeseries(row, validIndices); % Valid data
        invalidIndices = ~validIndices; % Indices to interpolate
        
        % Interpolate only for the invalid indices
        %interpolatedTimeseries(row, invalidIndices) = interp1(validTime, validData, t(invalidIndices), 'spline');
        interpolatedTimeseries(row, invalidIndices) = interp1(validTime, validData, t(invalidIndices), 'linear');
        
    end
    
    % Remove padded data
    timeseries = interpolatedTimeseries(:, paddingLength+1:end-paddingLength);
    %change motion mask to use all frames since you've already interpolated the data.
    motion_mask=ones(num_time_frames,1);
end

%% Start
for t = 1:num_time_frames
    disp(['Processing frame: ' num2str(t) '/' num2str(num_time_frames)]);
    
    if motion_mask(t) ==1
        if only_use_positive_data==1
            target_frame=double(timeseries(:,t)); % convert to double
            target_frame(timeseries(:,t)<0) =0;
        else
            target_frame = double(timeseries(:,t));
        end
        % solve the linear system to find contributions
        %x= eta_to_template_vox \ target_frame;
        %solve the non-negative least squares problem
        x=lsqnonneg(eta_to_template_vox,target_frame);
        
        %Calculate the modeled target frame based on  network contributions
        model_frame = eta_to_template_vox * x;
        
        %calculate the residul (unexplained proportion of the target frame)
        
        residual_vector = target_frame - model_frame;
        
        residual_norm = norm(residual_vector); % magnitude of the residual
        total_norm = norm (target_frame); %Magnitude of the target frame
        
        % calculate the proportion of each network and the residual
        
        
        if include_residuals==1
            network_proportions = (x / sum(x)) * (100 - (residual_norm / total_norm) *100);
            %residual_proportion = (residual_norm/total_norm *100);
            
        else
            %network_proportions = (norm(x) / total_norm) *100; % convert each weight to a percent
            network_proportions = (x / sum(x)) *100; % convert each weight to a percent
            
        end
        %Store the network proportions and residual proportions
        proportions_matrix(t,1:num_networks) = network_proportions;
        %fprintf('%s: %.2f%%\n', network_names{i}, network_proportions(i))
        if include_residuals ==1
            proportions_matrix(t,end) = (residual_norm / total_norm)*100;
        end
    end
end

%Get valid (non-censored) frames and time dor interpolation
valid_times = find(motion_mask);
valid_proportions = proportions_matrix(valid_times,:);

%use a cubis spline interpolation for missing frames
if do_spline_interpolation_on_nets ==1
    interpolated_proportions_matrix = proportions_matrix; % copy from the original matrix
    for i = 1:size(proportions_matrix,2) % loop over each network + residual;
        % do spline interpolation
        %interpolated_proportions_matrix(:,i) = spline(valid_times,valid_proportions(:,i),1:num_time_frames);
        %usepchip to avoid overshoots instead of spline.
        interpolated_proportions_matrix(:,i) = pchip(valid_times,valid_proportions(:,i),1:num_time_frames);
    end
end
%display the contribution percentage
fprintf('proportional contributions of each network template:\n');
% for i = 1:length(x)
%     %convert each network contribution to a proportion
%     proportion = (x(i) / sum(x) * 100 - residual_proportion); % scale by 100 %-residual
%     fprintf('%s: %.2f%%\n', network_names{i}, proportion)
% end

%display residual propotion
%fprintf ('Unexplained residual: %.2f%%\n', residual_proportion)

for t =1:image_limit
    disp(['frames = ' num2str(t)]);
    close all
%PLot data
time_vector = 1:num_time_frames; %replace with actual time vector if available;
figure();set(gcf,'Color','w');
if do_spline_interpolation_on_nets ==1
    a=area(time_vector, interpolated_proportions_matrix, 'Linestyle', 'none');
else
    a=area(time_vector, proportions_matrix, 'Linestyle', 'none');
    
end
xlabel('BOLD Frame');
ylabel('Network (%)')
title('Network Contributions Over Time');

% Customized legend with the network names and residual
if include_residuals ==1
    network_names_with_residual = [network_names, {'Unexplained Residual'}];
else
    network_names_with_residual=network_names;
end

legend(network_names_with_residual, 'Location','eastoutside');

%customize plot appearance
colormap(jet(num_networks +1)); % colormap
ylim([0, 100]);
netRGBs_no_missing = [
    255 0 0;
    0 0 153
    255 255 0
    0 255 0
    13 133 160
    50 50 50
    102 0 204
    102 255 255
    255 128 0
    178 102 153
    0 102 153
    102 255 102
    60 60 251
    200 200 200
    128 0 0
    ]/255;
%colororder(netRGBs);
for i =1:num_networks
    a(i).FaceColor = netRGBs_no_missing(i,:);
end
ylim([0, 100]);
xlim([t-50, t+50]); % 25 to make time "play by faster"
line([t t],[0 100],'Color','red','LineWidth',3)
%set(gca,'position',[0.025 0.11 0.95 0.8]);
%set(gcf,'Position',[33 343 1419 276]);
set(gcf,'Position',[33 343 973 200]);

print([output_mat_name num2str(t) '.png'], '-dpng', '-r300') % save image.
end



%total_time = size(timeseries,1)*TR;
for i= 1:size(proportions_matrix,1)
    max_value=max(proportions_matrix(i,:));
    main_net_columns(i)=find(proportions_matrix(i,:)==max_value);
end
% b=area(time_vector,main_net_columns);
% 
% for i =1:num_networks
%     b(i).FaceColor = netRGBs_no_missing(i,:);
% end

unique_nets = unique(main_net_columns);
n_states = length(unique_nets);

transition_counts = zeros(n_states,n_states);

%transitions include self-transitions
for t = 1:(length(main_net_columns) - 1)
    from_state = main_net_columns(t);
    to_state = main_net_columns(t+1);
    transition_counts(from_state, to_state) = transition_counts(from_state,to_state) + 1;
end
transition_matrix = transition_counts ./sum(transition_counts,2);
transition_matrix(isnan(transition_matrix)) =0; %handle nan rows
disp(transition_matrix)

figure(); set(gcf,'Color','w');
imagesc(transition_matrix);
xlabel('From Network');
ylabel('To Network');
colormap jet

figure(); set(gcf,'Color','w');
imagesc(main_net_columns)
colormap(netRGBs_no_missing)
xlim([0, 200]);

figure(); set(gcf,'Color','w');
[idx_kmeans,centroids] = kmeans(proportions_matrix,k);
imagesc(idx_kmeans'); colormap jet
xlim([0, 200]);

figure(); set(gcf,'Color','w');
total_time=sum(proportions_matrix);
proportion_time=sum(proportions_matrix)/sum(total_time);
%bar(categorical(network_names),proportion_time); % sorts alphabetically?
barx_names=categorical(network_names);
barx_names=reordercats(barx_names,network_names);
c = bar(barx_names,proportion_time); % sorts alphabetically?
c.FaceColor ='flat'; % Have to set this first before setting color.
for i=1:size(c.CData)
    c.CData(i,:) = netRGBs_no_missing(i,:);
end

%set(gca,'XTick',[]);

% calculate dwell time
% dwell_times = zeros(num_networks,1);
% for i =1:num_networks
%     current_state = unique_nets(i);
%     changes = [0, diff(main_net_columns ==current_state)];
%     segment_lengths = diff(find([changes,1] ==1));
%     dwell_times(i) = sum(segment_lengths);
% end
% 
% [a,b] = histc(main_net_columns,unique(main_net_columns));
% y = a(b);

save([output_mat_name '.mat'],'proportion_time','network_weights_filepath', 'dtseries_file_path',...
    'motion_mask','output_mat_name','only_use_positive_data',...
    'include_residuals','do_spline_interpolation_on_nets',...
    'do_spline_interpolation_on_dtseries',...
    'paddingLength','do_motion_censoring',...
    'use_cortex_only','transition_matrix');
disp('Done.')
end