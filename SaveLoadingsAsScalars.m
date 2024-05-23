function SaveLoadingsAsScalars(loadings_mat_file,template_dscalar_path,network_names_cells,wb_command,inherit_path,output_path)

%This function is designed to convert a .mat file as a dscalar file, one
%per network.

% arguements are:
% 1) a path to the .mat file from template matching
% 2) template dscalar file to save data into.
% 3) a list of network names, so that each dscalar has a unique name.
%    Format: 1 row where each cell in the row should contain the network name.
% 4) path to workbench command.

addpath(genpath('/home/faird/shared/code/external/utilities/gifti-1.6'))
addpath(genpath('/home/faird/shared/code/internal/utilities/Matlab_CIFTI'))
if exist('wb_command','var') ==1
else
    wb_command=settings.path_wb_c; %path to wb_command
end
try
    load(loadings_mat_file,'eta_to_template_vox','network_names','new_subject_labels');
    if exist('eta_to_template_vox','var') ==1
    else
        load(loadings_mat_file,'seed_matrix');
        eta_to_template_vox = seed_matrix;
    end
catch
    disp('eta_to_template_vox variable names not found. Trying to look for seematrix')
    try
        load(loadings_mat_file,'seed_matrix');
        eta_to_template_vox = seed_matrix;
    catch
        disp('unable to find 91282 x 16 data frame variable to use. exiting...');
        return
    end
    
end

if isempty(network_names_cells) ==1
    %use default network numbers
    disp('Using default network names loaded from .mat file.');
    %net_list = [1 2 3 5 7 8 9 10 11 12 13 14 15 16]; % hardcoded network assingments.
    network_names_cells =network_names;
else
    
end

template_cii = ciftiopen(template_dscalar_path,wb_command);
[matfilepath, matfilename]=fileparts(loadings_mat_file);
for i=1:length(network_names_cells)
    disp([network_names_cells{i} ' ' num2str(i)])
    if  i~=4 && i~=6 && i~=17
        this_scalar_data = eta_to_template_vox(:,i);
        template_cii.cdata=this_scalar_data;
        if inherit_path ==1
            ciftisave(template_cii,[matfilepath filesep matfilename '_' network_names_cells{i} '_loadings.dscalar.nii'],wb_command);
        else
            ciftisave(template_cii,[output_path filesep matfilename '_' network_names_cells{i} '_loadings.dscalar.nii'],wb_command);
            
        end
    else
    end
end
disp('Done making dscalars.')
end