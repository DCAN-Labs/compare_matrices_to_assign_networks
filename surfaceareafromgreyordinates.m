function [ network_surfarea, network_volume ] = surfaceareafromgreyordinates(Lmidthicknessfile,Rmidthicknessfile,dscalarwithassignments,outputname)

%This wrapper create a surface area dscalar

%% Add necessary paths
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

conc = strsplit(dscalarwithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments);
else
    dscalarwithassignments = {dscalarwithassignments};
end

%check to make sure that surface files exist
for i = 1:length(dscalarwithassignments)
    if exist(dscalarwithassignments{i}) == 0
        NOTE = ['Subject dscalar ' num2str(i) ' does not exist']
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

%check to make sure that surface files exist
for i = 1:length(Lmidthicknessfile)
    if exist(Lmidthicknessfile{i}) == 0
        NOTE = ['Subject L surface ' num2str(i) ' does not exist']
        disp(Lmidthicknessfile{i});
        return
    else
    end
end
disp('All Lsurface files exist continuing ...')

conc = strsplit(Rmidthicknessfile, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    Rmidthicknessfile = importdata(Rmidthicknessfile);
else
    Rmidthicknessfile = {Rmidthicknessfile};
end

%check to make sure that surface files exist
for i = 1:length(Rmidthicknessfile)
    if exist(Rmidthicknessfile{i}) == 0
        NOTE = ['Subject R surface ' num2str(i) ' does not exist']
        disp(Rmidthicknessfile{i});
        return
    else
    end
end
disp('All R surface files exist continuing ...')



newcii = ciftiopen(dscalarwithassignments{i},wb_command);
network_assignments = newcii.cdata;
num_networks = max(network_assignments);
surfscalar = ones(size(network_assignments,1),1);
surfscalar = surfscalar*8000;
newcii.cdata = surfscalar;
ciftisave(newcii, [support_folder '/' outputname '8000s_91282.dscalar.nii'],wb_command)

for i = 1:size(dscalarwithassignments) % loop through every subject
    newcii = ciftiopen(dscalarwithassignments{i},wb_command);
    network_assignments = newcii.cdata;
    
    cmd = ['wb_command -cifti-separate ' [support_folder '/' outputname '8000s_91282.dscalar.nii'] ' COLUMN -volume-all ' support_folder '/8000s_volume.nii'];
    system(cmd); clear cmd
    
    [~, A] = fileparts(Lmidthicknessfile{i});
    [~, B] = fileparts(A);
    C = [B '.vertexarea.func.gii'];
    cmd = ['wb_command -surface-vertex-areas '  Lmidthicknessfile{i} ' ' support_folder '/' C];
    disp(cmd)
    system(cmd); clear cmd
    
    [~, D] = fileparts(Rmidthicknessfile{i});
    [~, E] = fileparts(D);
    F = [E '.vertexarea.func.gii'];
    cmd = ['wb_command -surface-vertex-areas '  Rmidthicknessfile{i} ' ' support_folder '/' F];
    disp(cmd)
    system(cmd); clear cmd
    
    cmd = ['wb_command -cifti-create-dense-from-template ' [support_folder '/' outputname '8000s_91282.dscalar.nii'] ' ' [outputname 'surfaceareas.dscalar.nii'] ' -volume-all ' support_folder '/8000s_volume.nii -metric CORTEX_LEFT  ' support_folder '/' C ' -metric CORTEX_RIGHT ' support_folder '/' F];
    disp(cmd)
    system(cmd); clear cmd
    
    all_areas = ciftiopen([outputname 'surfaceareas.dscalar.nii'], wb_command);
    all_areas_vec = all_areas.cdata;
    
    
    for j = 1:num_networks
        net_indices = find(network_assignments == j);
        grey_SA = all_areas_vec(net_indices);
        sub_indices = find(grey_SA > 1000);
        surf_indices = find(grey_SA < 1000);
        network_volume(i,j) = size(sub_indices,1)*8;
        surf_SA = grey_SA(surf_indices);
        network_surfarea(i,j) = sum(surf_SA);
    end
end
disp('saving volumes and surface areas');
save( [pwd '/' outputname '.mat'], 'network_volume','network_surfarea')

disp('Done calculating surfaces areas and volumes');
end



