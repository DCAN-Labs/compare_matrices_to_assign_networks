function [all_areas_vec, network_surfarea, network_volume ] = surfaceareafromgreyordinates(Lmidthicknessfile_path,Rmidthicknessfile_path,output_only_greySA,dscalarwithassignments_path,outputname,output_folder,cleanupintermediatefiles,data_is_surface_only,restrict_to_ROI,ROI_vector)

%This wrapper create a surface area dscalar

% You'll need the subject's own midthickness files  and dscalar (dscalar.nii)of the network assignments.
% Lmidthicknessfile = a left mid-thickness file of the brain (.surf.gii) or (conc file of Lmidthckness files).
% Rmidthicknessfile = a right mid-thickness file of the brain (.surf.gii) or (conc file of Lmidthckness files).
% output_only_greySA = set to 1 to get the surface areas.  Set to 0 to debug.
% dscalarwithassignnments = vector of assignments for neural network numbers (1-18)   (dscalar.nii) or (conc file of dscalars)
% outputname = name of resulting file
% outputfolder = path to your output folder
% cleanupintermediatefiles = This script does a lot of file conversions, set to 1 to clean up. and 0 to inspect the intermediate files. 
% data_is_surface_only set to 1 if your data is surface only (i.e. 59412 grayordinates).  Set to 0 if you have a regular dscalar (i.e. 91282 grayordinates) 
% restrict_to_ROI = set to 1 if you want to restrict the surface area calculation to a specific ROI. set to 0 if you want to preform the surface area calculation across the whole brain.
% ROI_vector = optional. path to the an .mat file with a variable 'label_indices' vector where 1s and 0s where 1 is a grayordinate to keep and 0 is a grayordinate to ignore when doing the surface area calculation.  if you want to preform the surface area calculation across the whole brain, provide an empty string (e.g. '').

%example call using a conc file and and restricted ROI:
%surfaceareafromgreyordinates('/path/to/my_Lmidthickness.conc','/path/to/my_Rmidthickness.conc',1,'/path/to/my_dscalars.conc','my5subjectsoutputname','/path/to/myoutputfolder',1,0,1,'/path/to/myfrontal_cortex_indices.mat');
%example call for a single subject:
%surfaceareafromgreyordinates('/path/to/my_Lmidthickness.surf.gii','/path/to/my_Rmidthickness.surf.gii',1,'/path/to/my_network.dscalar.nii','myoutputname','/path/to/myoutputfolder',1,0,0,'');


%% Add necessary paths
this_code = which('surfaceareafromgreyordinates');
[code_dir,~] = fileparts(this_code);
[code_dir,~] = fileparts(code_dir);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

disp('Attempting to add neccesaary paths and functions.')
warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end

%add path to ft_read_cifti_mod so that you can get the left hemisphere  and
%right hemisphere mappings.
addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti')); % remove non-working gifti path included with MSCcodebase
warning('on')

wb_command=settings.path_wb_c; %path to wb_command
wb_command='/home/faird/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command';
%Lmidthicknessfile,Rmidthicknessfile,dscalarwithassignments

conc = strsplit(dscalarwithassignments_path, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments_path);
else
    dscalarwithassignments = {dscalarwithassignments_path};
end
tic
%check to make sure that surface files exist
% for i = 1:length(dscalarwithassignments)
%     if rem(i,100)==0
%         disp([' Validating file existence ' num2str(i)]);toc;
%     end
%     if exist(dscalarwithassignments{i}, 'file') == 0
%         disp(['Error Subject dscalar ' num2str(i) ' does not exist'])
%         disp(dscalarwithassignments{i});
%         return
%     else
%     end
% end
disp('All series files exist continuing ...')


conc = strsplit(Lmidthicknessfile_path, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    Lmidthicknessfile = importdata(Lmidthicknessfile_path);
else
    Lmidthicknessfile = {Lmidthicknessfile_path};
end
% check to make sure that surface files exist
% for i = 1:length(Lmidthicknessfile)
%     if rem(i,100)==0
%         disp([' Validating file existence ' num2str(i)]);toc;
%     end
%     if exist(Lmidthicknessfile{i}, 'file') == 0
%         disp(['Error Subject surface ' num2str(i) ' does not exist'])
%         disp(Lmidthicknessfile{i});
%         return
%     else
%     end
% end
disp('All series files exist continuing ...')


conc = strsplit(Rmidthicknessfile_path, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    Rmidthicknessfile = importdata(Rmidthicknessfile_path);
else
    Rmidthicknessfile = {Rmidthicknessfile_path};
end
% for i = 1:length(Rmidthicknessfile)
%     if rem(i,100)==0
%         disp([' Validating file existence ' num2str(i)]);toc;
%     end
%     if exist(Rmidthicknessfile{i}, 'file') == 0
%         disp(['Error Subject surface ' num2str(i) ' does not exist'])
%         disp(Rmidthicknessfile{i});
%         return
%     else
%     end
% end
disp('All series files exist continuing ...')

if isempty(ROI_vector)
else
    if ischar(ROI_vector) ==1
        load(ROI_vector);
        ROI_vector=label_indices;
    end
    
end
%Make a blank dscalar with 8000 as all the values using the 1st cifti.
newcii = ciftiopen(dscalarwithassignments{1},wb_command);
network_assignments = newcii.cdata;

newcii_ft = ft_read_cifti_mod(dscalarwithassignments{1});
cii_struct_labels_indices = find(newcii_ft.brainstructure>0);
cii_struct_labels = newcii_ft.brainstructure(cii_struct_labels_indices);
                
indices_left_log = cii_struct_labels==1;
indices_right_log = cii_struct_labels==2;
                
num_networks = max(network_assignments); %check that subject 1 has all networks.
surfscalar = ones(size(network_assignments,1),1);
surfscalar = surfscalar*8000;
newcii.cdata = surfscalar;
if data_is_surface_only ==1
    ciftisave(newcii, [output_folder filesep outputname '_59412.dscalar.nii'],wb_command)
    
else
    ciftisave(newcii, [output_folder filesep outputname '8000s_91282.dscalar.nii'],wb_command)
    
end

%preallocate
network_surfarea = zeros(size(dscalarwithassignments,1),num_networks);
network_volume = zeros(size(dscalarwithassignments,1),num_networks);
network_L_hemi_surfarea = zeros(size(dscalarwithassignments,1),num_networks);
network_R_hemi_surfarea = zeros(size(dscalarwithassignments,1),num_networks);

for i = 1:size(dscalarwithassignments) % loop through every subject
    disp(['Getting surface area and volume for subject: ' num2str(i)])
    
    newcii = ciftiopen(dscalarwithassignments{i},wb_command);
    network_assignments = newcii.cdata;
    
    if data_is_surface_only
    else
        cmd = [wb_command ' -cifti-separate ' [output_folder filesep outputname '8000s_91282.dscalar.nii'] ' COLUMN -volume-all ' output_folder filesep '8000s_volume.nii'];
        system(cmd); clear cmd
    end
    
    [~, A] = fileparts(Lmidthicknessfile{i});
    [~, B] = fileparts(A);
    C = [B '.vertexarea.func.gii'];
    cmd = [wb_command ' -surface-vertex-areas '  Lmidthicknessfile{i} ' ' output_folder filesep C];
    disp(cmd)
    system(cmd);
    clear cmd
    
    [~, D] = fileparts(Rmidthicknessfile{i});
    [~, E] = fileparts(D);
    F = [E '.vertexarea.func.gii'];
    cmd = [wb_command ' -surface-vertex-areas '  Rmidthicknessfile{i} ' ' output_folder filesep F];
    disp(cmd)
    system(cmd);
    clear cmd
    
    if data_is_surface_only ==1
        cmd = [wb_command ' -cifti-create-dense-from-template ' [output_folder filesep outputname '_59412.dscalar.nii'] ' ' [outputname 'surfaceareas.dscalar.nii'] ' -metric CORTEX_LEFT  ' output_folder filesep C ' -metric CORTEX_RIGHT ' output_folder filesep F];
        disp(cmd)
        system(cmd); clear cmd
    else
        cmd = [wb_command ' -cifti-create-dense-from-template ' [output_folder filesep outputname '8000s_91282.dscalar.nii'] ' ' [outputname 'surfaceareas.dscalar.nii'] ' -volume-all ' output_folder filesep '8000s_volume.nii -metric CORTEX_LEFT  ' output_folder filesep C ' -metric CORTEX_RIGHT ' output_folder filesep F];
        disp(cmd)
        system(cmd); clear cmd
    end
    
    all_areas = ciftiopen([outputname 'surfaceareas.dscalar.nii'], wb_command);
    all_areas_vec = all_areas.cdata;
        all_areas_vec_orig = all_areas.cdata;

    if output_only_greySA   ==1
        network_assignment_filetype = strsplit(dscalarwithassignments{i}, '.');
        cifti_type = char(network_assignment_filetype(end-1));
        if strcmp('dtseries',cifti_type) == 1
            overlap =1;
        else
            overlap =0;
        end

        
        if overlap ==1
            for j = 1:size(network_assignments,2)
                %for j = 1:num_networks
                net_indices = find(network_assignments(:,j) ~= 0);
                net_indices_log = network_assignments == j;
                
                net_and_left = net_indices_log & indices_left_log; %logical sum of hemisphere and network
                net_indices_left = find(net_and_left==1);
                net_and_right = net_indices_log & indices_right_log; %logical sum of hemisphere and network
                net_indices_right = find(net_and_right==1);
                
                grey_SA = all_areas_vec(net_indices);
                sub_indices = find(grey_SA > 1000);
                surf_indices = find(grey_SA < 1000);
                if data_is_surface_only ==1
                    network_volume(i,j) = 0;
                else                    
                    network_volume(i,j) = size(sub_indices,1)*8;
                end
                total_net_SA = grey_SA(surf_indices);
                network_surfarea(i,j) = sum(total_net_SA);
                
                Lhemi_net_SA = all_areas_vec(net_indices_left);
                Rhemi_net_SA = all_areas_vec(net_indices_right);
                network_L_hemi_surfarea(i,j) = sum(Lhemi_net_SA);
                network_R_hemi_surfarea(i,j) = sum(Rhemi_net_SA);
                
            end
            
        else % run single network assingment
            
            for j = 1:num_networks
                if restrict_to_ROI ==0
                    net_indices = find(network_assignments == j);
                    net_indices_log = network_assignments == j;
                    net_and_left = net_indices_log & indices_left_log; %logical sum of hemisphere and network
                    net_and_right = net_indices_log & indices_right_log; %logical sum of hemisphere and network
                    
                else
                    net_indices = find(network_assignments((ROI_vector)) == j);
                    net_indices_log = network_assignments((ROI_vector)) == j;
                    net_and_left = net_indices_log & indices_left_log(ROI_vector); %logical sum of hemisphere and network
                    net_and_right = net_indices_log & indices_right_log(ROI_vector); %logical sum of hemisphere and network
                    
                end
                
                net_indices_left = find(net_and_left==1);
                net_indices_right = find(net_and_right==1);
                if restrict_to_ROI ==1
                    all_areas_vec =all_areas_vec_orig(ROI_vector);
                end
                grey_SA = all_areas_vec(net_indices);
                
                sub_indices = find(grey_SA > 1000);
                surf_indices = find(grey_SA < 1000); %if surface only data is incorrectly indicated, then this function should still work, because the surface area of the addional grayordinates is zero.
                if data_is_surface_only ==1
                    network_volume(i,j) = 0;
                else
                    network_volume(i,j) = size(sub_indices,1)*8;
                end
                
                total_net_SA = grey_SA(surf_indices);
                Lhemi_net_SA = all_areas_vec(net_indices_left);
                Rhemi_net_SA = all_areas_vec(net_indices_right);
                
                %if restrict_to_ROI ==1
                %network_surfarea_restricted(i,j) = sum(total_net_SA(ROI_vector));
                %end
                network_surfarea(i,j) = sum(total_net_SA);
                
                network_L_hemi_surfarea(i,j) = sum(Lhemi_net_SA);
                network_R_hemi_surfarea(i,j) = sum(Rhemi_net_SA);
            end
        end
    else
    end
    
    if cleanupintermediatefiles ==1
        cmd = ['rm ' output_folder filesep C ' ' output_folder filesep F ' ' ];
        system(cmd);
    else
    end
    
    
end

if cleanupintermediatefiles ==1
    cmd = ['rm ' output_folder filesep outputname '8000s_91282.dscalar.nii ' output_folder filesep '8000s_volume.nii'];
    system(cmd)
else
end


for i= 1:size(network_surfarea,1)
network_proportion(i,:) = network_surfarea(i,:)/(sum(network_surfarea(i,:)));
end


if output_only_greySA  ==1
    disp('saving volumes and surface areas');
    if restrict_to_ROI ==1
        save( [output_folder filesep outputname '.mat'], 'network_volume','network_surfarea','network_L_hemi_surfarea','network_R_hemi_surfarea','dscalarwithassignments','Lmidthicknessfile','Rmidthicknessfile','network_proportion','restrict_to_ROI','ROI_vector')
        
    else
        save( [output_folder filesep outputname '.mat'], 'network_volume','network_surfarea','network_L_hemi_surfarea','network_R_hemi_surfarea','dscalarwithassignments','Lmidthicknessfile','Rmidthicknessfile','network_proportion')
    end
else
end

disp('Done calculating surfaces areas and volumes');
end



