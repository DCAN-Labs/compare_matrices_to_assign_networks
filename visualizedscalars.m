function visualizedscalars(dscalarswithassignments,outputname,output_map_type, plot_results)

%This code loads in a conc of dcalars to visualize them for all subjects.
% It also can calculates the probability of a network assingment from the
% list of subjects that you've provided.

% dscalar_with_assingment = a .conc file of all your subject's dscalars (1 path per row).
% Save  percentages = will save a file (dscalar for each network).

%Visualization Paremeters:
%Downsample_scalar = if true, the dscalar with be down sampled (sampled
%according to the DS_factor.  If false, it will not downsample data for visualization.
%DS_factor = 50; %downsample factor.  Reduce the 91282 vector by this factor to reduce the load on matlab visualiztion tools.
%(e.g. DS = 2, sample every other greyordinate.  visualization will have
%45641 data points per subject).

DS_factor = 50; %downsample factor.  Reduce the 91282 vector by this factor to reduce the load on matlab visualiztion tools.
%(e.g. DS = 2, sample every other greyordinate.  visualizetion will have
%45641 data points per subject).
downsample_scalar = 1;
save_results = 1;
check_for_nets_greater_than_16 = 0; % Set to 0 if running infomap. Otherwise the code will warn you that for each subject  than 16.
try_resort_on_subject = 0;
if isnumeric(DS_factor)==1
else
    DS_factor = str2num(DS_factor);
end

if isnumeric(save_results)==1
else
    save_results = str2num(save_results);
end

switch output_map_type
    
    case 'calc_percentage'
        
        calc_percentage =1;
        number_of_networks =0;
        calc_probability = 0;
        
    case 'number_of_networks'
        
        calc_percentage =0;
        number_of_networks=1;
        calc_probability = 0;
        
    case 'calc_probability'
        
        calc_percentage =0;
        number_of_networks =0;
        calc_probability = 1;

     case 'calc_mode'
        
        calc_percentage =0;
        number_of_networks =0;
        calc_probability = 1;
        calc_mode =1;
        
    otherwise
        disp('output map type not supported. check your inputs.')
end
%load colormap
load('/panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/PowerColorMap.mat')
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};
conc = strsplit(dscalarswithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarswithassignments = importdata(dscalarswithassignments);
else
    dscalarswithassignments = {dscalarswithassignments};
end

num_orig_files = length(dscalarswithassignments);
network_assignment_filetype = strsplit(dscalarswithassignments{1}, '.');
cifti_type = char(network_assignment_filetype(end-1));
if strcmp('dtseries',cifti_type) == 1
    overlap =1;
else
    overlap =0;
end
tic

%check to make sure that surface files exist
found_files_total=0; missing_files_total =0;% make a "found files counter" 
for i = 1:length(dscalarswithassignments)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(dscalarswithassignments{i}, 'file') == 0
        disp(['Error Subject dscalar ' num2str(i) ' does not exist'])
        disp(dscalarswithassignments{i});
        missing_files_total = missing_files_total+1;
        missing_files_indx(missing_files_total) = i;
        missing_files{missing_files_total} = dscalarswithassignments{i};
        %return
    else
        found_files_total = found_files_total+1;
        found_files_indx(found_files_total) = i;
        found_files{found_files_total} = dscalarswithassignments{i};
    end
end

if found_files_total ==length(dscalarswithassignments)
    disp('All series files exist continuing ...')
else
    disp('WARNING: Not all files were found.')
    disp(['Expected to find: ' num2str(length(dscalarswithassignments))])
    disp(['Files found: ' num2str(found_files_total)])
    disp(['Files missing: ' num2str(missing_files_total)])
    prompt = 'Continue with only found files? [Y/N]: ';
    str = input(prompt, 's');
    if strcmp(str,'y')==1 || strcmp(str,'Y')==1 || strcmp(str,'yes')==1 || strcmp(str,'YES')==1 || strcmp(str,'Yes')==1 || isempty(str)
        disp('Using only found files.')
        dscalarswithassignments = found_files;
    else
        return
    end
end

%% Adding paths for this function
this_code = which('visualizedscalars');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
%support_folder=[pwd '/support_files'];
if ~isdeployed
    addpath(genpath(support_folder));
end
settings=settings_comparematrices;%
%settings.path';
np=size(settings.path,2);

if ~isdeployed
    warning('off') %supress addpath warnings to nonfolders.
    for i=1:np
        addpath(genpath(settings.path{i}));
    end
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
wb_command=settings.path_wb_c; %path to wb_command
warning('on')



%use the first 
disp('Using the first file to infer the number of greyordinates.')
scalar_temp = ciftiopen(dscalarswithassignments{1},wb_command);
grey_size=size(scalar_temp.cdata,1);
switch grey_size
    
    case 91282
        disp('Number of greyordinates is 91282.')
        surface_only =0;
    case 59412
        disp('Number of greyordinates is 59412.')
        surface_only =1;
    otherwise
        disp('Number of greyordinates is neither 91282 nor 59412.  You may run into saving issues...')    
        surface_only =0;
end

if (overlap ==1 && calc_percentage == 1) || (overlap ==1 && calc_probability == 1)
    scalar_array = zeros(grey_size,length(dscalarswithassignments),length(network_names)); %preallocate
else
    scalar_array = zeros(grey_size,length(dscalarswithassignments)); %preallocate
end

for i=1:length(dscalarswithassignments)
%for i=1:5 % for debugging
    disp(i)
    scalar_temp = ciftiopen(dscalarswithassignments{i},wb_command);
    scalar=scalar_temp.cdata;
    if overlap == 1
        if calc_percentage ==1
            scalar_array(:,i,:) =scalar; %use scalar name despite it being a matrix.
        elseif calc_probability == 1
            scalar_array(:,i,:) =scalar; %use scalar name despite it being a matrix.
        elseif  number_of_networks ==1
            islabeled = scalar ~=0;
            scalar_array(:,i) = sum(islabeled,2);
        elseif calc_mode ==1
           scalar_array(:,i,:) =scalar; %use scalar name despite it being a matrix.
        end
    else
        scalar_array(:,i) = scalar;
    end
    
end
if check_for_nets_greater_than_16 ==1
%check for values greater than 16.
for i=1:length(dscalarswithassignments)
    if (overlap ==1 && calc_percentage == 1) || (overlap ==1 && calc_probability == 1)
        for j=1:length(network_names)
            isgreaterthan16 = scalar_array(:,i,j) > 16;
            if sum(isgreaterthan16) ~= 0
                disp([dscalarswithassignments{i}])
            else
            end
        end
    else
        isgreaterthan16 = scalar_array(:,i) > 16;
        if sum(isgreaterthan16) ~= 0
            disp([dscalarswithassignments{i}])
        else
        end
    end
end
end
%disp('extract number of unique networks from subject 1.')
%num_networks=unique(scalar_array(:,1));

%open a file to write for saving
if surface_only ==1
    temp_file=ciftiopen('/panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii',wb_command);
else
    temp_file=ciftiopen('/panfs/roc/groups/8/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii',wb_command);
end
length(network_names)
size(scalar_array,1)
if overlap == 1
    if calc_percentage == 1 || calc_probability == 1
        for i=1:length(network_names) % i is network names
            disp(i)
            if  i~=4 && i~=6
                allthisnets=squeeze(scalar_array(:,:,i));
                %                 for j=1:size(scalar_array,1) % j is  voxel
                %                     if rem(j,5000)==0
                %                         %disp([' Calculating voxel proportion ' num2str(j)]);
                %                     end
                %                     is_in_network = find(scalar_array(j,:,i) == i);
                is_in_network = allthisnets == i;
                num_subjects_network_true = sum(is_in_network,2);
                if calc_percentage == 1
                    network_percentage(:,i) = (num_subjects_network_true/length(dscalarswithassignments))*100;
                else
                    network_percentage(:,i) = (num_subjects_network_true/length(dscalarswithassignments));
                end
                %                 end
                if save_results ==1
                    temp_file.cdata=network_percentage(:,i);
                    disp('Saving percentages/probabilty for each network.')
                    if calc_percentage == 1
                        ciftisave(temp_file,[outputname '_' network_names{i} '_network_percentage.dscalar.nii'],wb_command);
                    else
                        ciftisave(temp_file,[outputname '_' network_names{i} '_network_probability.dscalar.nii'],wb_command);
                    end
                    if found_files_total ==num_orig_files
                    else
                        disp('Saving list of found files.')
                        T = table(found_files,'VariableNames',{'found_files'});
                        writetable(T,[outputname '_found_files.txt']);
                    end
                end
            end
        end
    else %calculate number of networks
        avg_num_networks =   mean(scalar_array,2);
        if save_results ==1
            temp_file.cdata=  avg_num_networks;
            disp('Saving average number of networks.')
            whole_brain_number_of_nets = mean(scalar_array,1);
            whole_brain_number_of_nets = whole_brain_number_of_nets'; 
            integration_zonesscalarpath = '/panfs/roc/groups/3/rando149/shared/projects/ABCD_net_template_matching/ABCD_number_of_nets/ABCD_GRP1_overlap_number_of_nets_avg_number_of_network_2.2_thres_sz60_clusters.dscalar.nii';
            intcii = ciftiopen(integration_zonesscalarpath,wb_command);
            intmask = intcii.cdata;
            intmask = intmask>0; %because this was built from a label file, binarize to make a mask.
            intscalar=scalar_array(intmask,:);
            integration_zone_number_of_nets=mean(intscalar,1);
            integration_zone_number_of_nets = integration_zone_number_of_nets';
            save([outputname '.mat'],'whole_brain_number_of_nets','integration_zone_number_of_nets');
            ciftisave(temp_file,[outputname '_avg_number_of_networks.dscalar.nii'],wb_command);
        end
    end
else % files are dscalars.
    for i=1:length(network_names)
        disp(i)
        if  i~=4 && i~=6
            for j=1:size(scalar_array,1)
                if rem(j,5000)==0
                    %disp([' Calculating voxel proportion ' num2str(j)]);
                end
                is_in_network = find(scalar_array(j,:) == i);
                num_subjects_network_true = length(is_in_network);
                if calc_percentage ==1
                    network_percentage(i,j) = (num_subjects_network_true/length(dscalarswithassignments))*100;
                else
                    network_percentage(i,j) = (num_subjects_network_true/length(dscalarswithassignments));
                end
            end
            if save_results ==1
                temp_file.cdata=network_percentage(i,:)';
                disp('Saving percentages for each network.')
                if calc_percentage == 1
                    
                    ciftisave(temp_file,[outputname '_' network_names{i} '_network_percentage.dscalar.nii'],wb_command);
                else
                    ciftisave(temp_file,[outputname '_' network_names{i} '_network_probability.dscalar.nii'],wb_command);
                end
                if found_files_total ==num_orig_files
                else
                    disp('Saving list of found files.')
                    T = table(found_files,'VariableNames',{'found_files'});
                    writetable(T,[outputname '_found_files.txt']);
                end
            end
        end
    end
end
large_scalar_array = scalar_array;

if plot_results ==1
    if downsample_scalar ==1
        scalar_array = scalar_array(1:DS_factor:end,:);
        imagesc(scalar_array)
        if overlap == 0
            
            colormap(mymap)
        end
        
    else
        imagesc(scalar_array)
        if overlap == 0
            colormap(mymap)
        end
    end
    
    
    for i=1:length(dscalarswithassignments)
        sorted_array(:,i) = sort(scalar_array(:,i));
    end
    
    figure()
    imagesc(sorted_array)
    if overlap == 0
        colormap(mymap)
    end
    %try this for sorting
    %[~,sorted_single] = sort(scalar_array(:,1));
    if overlap == 0
        if try_resort_on_subject ==1
            [~,sorted_single] = sort(temp_file.cdata(:,1));
            for i=1:length(dscalarswithassignments)
                one_subject_sorted_array(:,i) = scalar_array(sorted_single,i);
            end
            figure()
            imagesc(one_subject_sorted_array)
            if overlap == 0
                colormap(mymap)
            end
            
            for i=1:length(dscalarswithassignments)
                num_default = find(scalar_array(:,i) == 1);
                all_num_default(i) = size(num_default,1);
            end
            rearranged_sorted_array = sorted_array(:,all_num_default);
            [~, idxb] = sort(all_num_default);
            rearranged_sorted_array = sorted_array(:,idxb);
            
            figure()
            imagesc(rearranged_sorted_array)
            if overlap == 0
                colormap(mymap)
            end
        else
        end
    end
end
disp('Done making network dscalars.')
