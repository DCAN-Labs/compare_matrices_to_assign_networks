function visualizedscalars(dscalarwithassignments)


%This code loads in a conc of dcalars to visualize them for all subjects.

DS_factor = 50; %downsample factor.  Reduce the 91282 vector by this factor to reduce the load on matlab visualiztion tools.
%(e.g. DS = 2, sample every other greyordinate.  visualizetion will have
%45641 data points per subject).
downsample_scalar =1;
save_percentages =1;
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};


conc = strsplit(dscalarwithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments);
else
    dscalarwithassignments = {dscalarwithassignments};
end

%check to make sure that surface files exist
for i = 1:length(dscalarwithassignments)
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
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
wb_command=settings.path_wb_c; %path to wb_command
warning('on')

scalar_array = zeros(91282,length(dscalarwithassignments)); %preallocate

for i=1:length(dscalarwithassignments)
    
    scalar_temp = ciftiopen(dscalarwithassignments{i},wb_command);
    scalar=scalar_temp.cdata;
    scalar_array(:,i) = scalar;
    
end

%check for values greater than 16.
for i=1:length(dscalarwithassignments)
    isgreaterthan16 = scalar_array(:,i) > 16;
    if sum(isgreaterthan16) ~= 0
    disp([dscalarwithassignments{i}])
    else
    end
end

disp('extract number of unique networks from subject 1.')
num_networks=unique(scalar_array(:,1));

%open a file to write for saving
    temp_file=ciftiopen('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/91282_Greyordinates.dscalar.nii',wb_command);
    for i=1:length(network_names)
        if i==4 || i==6
            continue
        end
        for j=1:size(scalar_array,1)
            if rem(j,5000)==0
                disp([' Calculating voxel proportion ' num2str(j)]);
            end
            is_in_network = find(scalar_array(j,:) == i);
            num_subjects_network_true = length(is_in_network);
            network_percentage(i,j) = (num_subjects_network_true/length(dscalarwithassignments))*100;
            
        end
        if save_percentages ==1
            temp_file.cdata=network_percentage(i,:)';
            disp('Saving percentages for each network.')
            ciftisave(temp_file,[network_names{i} '_network_percentage.dscalar.nii'],wb_command);
        else
        end
    end
    
large_scalar_array = scalar_array;  

%plot results
load('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/support_files/PowerColorMap.mat')


if downsample_scalar ==1
scalar_array = scalar_array(1:DS_factor:end,:);
imagesc(scalar_array)
colormap(mymap)


else
    imagesc(scalar_array)
    colormap(mymap)
end

    
for i=1:length(dscalarwithassignments) 
    sorted_array(:,i) = sort(scalar_array(:,i)); 
end
figure()
imagesc(sorted_array)
    colormap(mymap)
    
%try this for sorting
%[~,sorted_single] = sort(scalar_array(:,1));
[~,sorted_single] = sort(temp_file.cdata(:,1));
for i=1:length(dscalarwithassignments) 
    one_subject_sorted_array(:,i) = scalar_array(sorted_single,i); 
end
figure()
imagesc(one_subject_sorted_array)
    colormap(mymap)

    
for i=1:length(dscalarwithassignments)
    num_default = find(scalar_array(:,i) == 1);
    all_num_default(i) = size(num_default,1);
end
    rearranged_sorted_array = sorted_array(:,all_num_default);
    [~, idxb] = sort(all_num_default);
    rearranged_sorted_array = sorted_array(:,idxb);
  
    figure()
    imagesc(rearranged_sorted_array)
    colormap(mymap)

