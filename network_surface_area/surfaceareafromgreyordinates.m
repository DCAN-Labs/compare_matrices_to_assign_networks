function [all_areas_vec, network_surfarea, network_volume ] = surfaceareafromgreyordinates(Lmidthicknessfile,Rmidthicknessfile,output_only_greySA,dscalarwithassignments,outputname,output_folder,cleanupintermediatefiles,data_is_surface_only)

%This wrapper create a surface area dscalar

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

%Lmidthicknessfile,Rmidthicknessfile,dscalarwithassignments

conc = strsplit(dscalarwithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments);
else
    dscalarwithassignments = {dscalarwithassignments};
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


conc = strsplit(Lmidthicknessfile, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    Lmidthicknessfile = importdata(Lmidthicknessfile);
else
    Lmidthicknessfile = {Lmidthicknessfile};
end
% check to make sure that surface files exist
for i = 1:length(Lmidthicknessfile)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(Lmidthicknessfile{i}, 'file') == 0
        disp(['Error Subject surface ' num2str(i) ' does not exist'])
        disp(Lmidthicknessfile{i});
        return
    else
    end
end
disp('All series files exist continuing ...')


conc = strsplit(Rmidthicknessfile, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    Rmidthicknessfile = importdata(Rmidthicknessfile);
else
    Rmidthicknessfile = {Rmidthicknessfile};
end
for i = 1:length(Rmidthicknessfile)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(Rmidthicknessfile{i}, 'file') == 0
        disp(['Error Subject surface ' num2str(i) ' does not exist'])
        disp(Rmidthicknessfile{i});
        return
    else
    end
end
disp('All series files exist continuing ...')


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
    ciftisave(newcii, [output_folder '/' outputname '_59412.dscalar.nii'],wb_command)
    
else
    ciftisave(newcii, [output_folder '/' outputname '8000s_91282.dscalar.nii'],wb_command)
    
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
        cmd = ['wb_command -cifti-separate ' [output_folder '/' outputname '8000s_91282.dscalar.nii'] ' COLUMN -volume-all ' output_folder '/8000s_volume.nii'];
        system(cmd); clear cmd
    end
    
    [~, A] = fileparts(Lmidthicknessfile{i});
    [~, B] = fileparts(A);
    C = [B '.vertexarea.func.gii'];
    cmd = ['wb_command -surface-vertex-areas '  Lmidthicknessfile{i} ' ' output_folder '/' C];
    disp(cmd)
    system(cmd);
    clear cmd
    
    [~, D] = fileparts(Rmidthicknessfile{i});
    [~, E] = fileparts(D);
    F = [E '.vertexarea.func.gii'];
    cmd = ['wb_command -surface-vertex-areas '  Rmidthicknessfile{i} ' ' output_folder '/' F];
    disp(cmd)
    system(cmd);
    clear cmd
    
    if data_is_surface_only ==1
        cmd = ['wb_command -cifti-create-dense-from-template ' [output_folder '/' outputname '_59412.dscalar.nii'] ' ' [outputname 'surfaceareas.dscalar.nii'] ' -metric CORTEX_LEFT  ' output_folder '/' C ' -metric CORTEX_RIGHT ' output_folder '/' F];
        disp(cmd)
        system(cmd); clear cmd
    else
        cmd = ['wb_command -cifti-create-dense-from-template ' [output_folder '/' outputname '8000s_91282.dscalar.nii'] ' ' [outputname 'surfaceareas.dscalar.nii'] ' -volume-all ' output_folder '/8000s_volume.nii -metric CORTEX_LEFT  ' output_folder '/' C ' -metric CORTEX_RIGHT ' output_folder '/' F];
        disp(cmd)
        system(cmd); clear cmd
    end
    
    all_areas = ciftiopen([outputname 'surfaceareas.dscalar.nii'], wb_command);
    all_areas_vec = all_areas.cdata;
    
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
                net_indices = find(network_assignments == j);
                net_indices_log = network_assignments == j;
                
                net_and_left = net_indices_log & indices_left_log; %logical sum of hemisphere and network
                net_indices_left = find(net_and_left==1);
                net_and_right = net_indices_log & indices_right_log; %logical sum of hemisphere and network
                net_indices_right = find(net_and_right==1);  
                
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
                
                network_surfarea(i,j) = sum(total_net_SA);
                network_L_hemi_surfarea(i,j) = sum(Lhemi_net_SA);
                network_R_hemi_surfarea(i,j) = sum(Rhemi_net_SA);
            end
        end
    else
    end
    
    if cleanupintermediatefiles ==1
        cmd = ['rm ' output_folder '/' C ' ' output_folder '/' F ' ' ];
        system(cmd);
    else
    end
    
    
end

if cleanupintermediatefiles ==1
    cmd = ['rm ' output_folder '/' outputname '8000s_91282.dscalar.nii ' output_folder '/8000s_volume.nii'];
    system(cmd)
else
end

if output_only_greySA  ==1
    disp('saving volumes and surface areas');
    save( [output_folder '/' outputname '.mat'], 'network_volume','network_surfarea','network_L_hemi_surfarea','network_R_hemi_surfarea')
else
end

disp('Done calculating surfaces areas and volumes');
end



