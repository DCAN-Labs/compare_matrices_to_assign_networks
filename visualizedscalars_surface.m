function visualizedscalars_surface(dscalarwithassignments,outputname,output_map_type)

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

DS_factor = 50; %downsample factor.  Reduce the 59412 vector by this factor to reduce the load on matlab visualiztion tools.
%(e.g. DS = 2, sample every other greyordinate.  visualizetion will have
%45641 data points per subject).
downsample_scalar =1;
save_percentages =1;

if isnumeric(DS_factor)==1
else
    DS_factor = str2num(DS_factor);
end

if isnumeric(save_percentages)==1
else
    save_percentages = str2num(save_percentages);
end

switch output_map_type
    
    case 'calc_percentage'
        
        calc_percentage =1;
        number_of_networks =0;
    case 'number_of_networks'
        
        calc_percentage =0;
        number_of_networks=1;
        
    otherwise
        disp('output map type not supported. check your inputs.')
end
%load colormap
load('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/PowerColorMap.mat')
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};
conc = strsplit(dscalarwithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments);
else
    dscalarwithassignments = {dscalarwithassignments};
end

network_assignment_filetype = strsplit(dscalarwithassignments{1}, '.');
cifti_type = char(network_assignment_filetype(end-1));
if strcmp('dtseries',cifti_type) == 1
    overlap =1;
else
    overlap =0;
end
tic
%check to make sure that surface files exist
for i = 1:length(dscalarwithassignments)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(dscalarwithassignments{i}, 'file') == 0
        disp(['Error Subject dscalar ' num2str(i) ' does not exist'])
        disp(dscalarwithassignments{i});
        return
    else
    end
end
disp('All series files exist continuing ...')


%% Adding paths for this function
this_code = which('visualizedscalars');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
%support_folder=[pwd '/support_files'];
if ~isdeployed
    addpath(genpath(support_folder));
end
settings=settings_comparematrices;%
settings.path';
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

if overlap ==1 && calc_percentage == 1
    scalar_array = zeros(59412,length(dscalarwithassignments),length(network_names)); %preallocate
else
    scalar_array = zeros(59412,length(dscalarwithassignments)); %preallocate
end

for i=1:length(dscalarwithassignments)
    disp(i)
    scalar_temp = ciftiopen(dscalarwithassignments{i},wb_command);
    scalar=scalar_temp.cdata;
    if overlap == 1
        if calc_percentage ==1
            scalar_array(:,i,:) =scalar; %use scalar name despite it being a matrix.
        else
            islabeled = scalar ~=0;
            scalar_array(:,i) = sum(islabeled,2);
        end
    else
        scalar_array(:,i) = scalar;
    end
    
end

%check for values greater than 16.
for i=1:length(dscalarwithassignments)
    if overlap ==1 && calc_percentage == 1
        for j=1:length(network_names)
            isgreaterthan16 = scalar_array(:,i,j) > 16;
            if sum(isgreaterthan16) ~= 0
                disp([dscalarwithassignments{i}])
            else
            end
        end
    else
        isgreaterthan16 = scalar_array(:,i) > 16;
        if sum(isgreaterthan16) ~= 0
            disp([dscalarwithassignments{i}])
        else
        end
    end
end

disp('extract number of unique networks from subject 1.')
%num_networks=unique(scalar_array(:,1));

%open a file to write for saving
temp_file=ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/91282_Greyordinates_surf_only.dscalar.nii',wb_command);
length(network_names)
size(scalar_array,1)
if overlap == 1
    if calc_percentage == 1
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
                     network_percentage(:,i) = (num_subjects_network_true/length(dscalarwithassignments))*100;
%                     
%                 end
                if save_percentages ==1
                    temp_file.cdata=network_percentage(:,i);
                    disp('Saving percentages for each network.')
                    ciftisave(temp_file,[outputname '_' network_names{i} '_network_percentage.dscalar.nii'],wb_command);
                end
            end
        end
    else %calculate number of networks
        avg_num_networks =   mean(scalar_array,2);
        if save_percentages ==1
            temp_file.cdata=  avg_num_networks;
            disp('Saving average number of networks.')
            ciftisave(temp_file,[outputname '_avg_number_of_network.dscalar.nii'],wb_command);
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
                network_percentage(i,j) = (num_subjects_network_true/length(dscalarwithassignments))*100;
                
            end
            if save_percentages ==1
                temp_file.cdata=network_percentage(i,:)';
                disp('Saving percentages for each network.')
                ciftisave(temp_file,[outputname '_' network_names{i} '_network_percentage.dscalar.nii'],wb_command);
            end
        end
    end
end
large_scalar_array = scalar_array;

%plot results


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


for i=1:length(dscalarwithassignments)
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
    [~,sorted_single] = sort(temp_file.cdata(:,1));
    for i=1:length(dscalarwithassignments)
        one_subject_sorted_array(:,i) = scalar_array(sorted_single,i);
    end
    figure()
    imagesc(one_subject_sorted_array)
    if overlap == 0
        colormap(mymap)
    end
    
    for i=1:length(dscalarwithassignments)
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
end
disp('Done.')
