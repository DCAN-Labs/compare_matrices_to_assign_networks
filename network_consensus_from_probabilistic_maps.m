function network_consensus_from_probabilistic_maps(path_to_probability_maps_list,run_locally,output_name,output_folder,run_clean_up_script)

%This function loads a list of probabilistic maps and generates a
%whole-brain parecelation schema based on the highest region with the
%highest probability.

%inputs are: 
%path_to_probability_maps_list= a list of dcscalars that are probabilistic
%maps.
%assumed network order in conc file. =  {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};
%run_locally, set to 1 if your Robert and running this on your local
%computer. Otherwise set to 0. (This option points to cifti dependencies.)
%output_name= full output name. i.e. /path/to/my/dscalar.nii
%run_clean_up_script= remove small islands of networks smaller than 20
%grayordinates.

%% Adding paths for this function
this_code = which('brain_gerrymander');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir filesep 'support_folder']; %find support files in the code directory.
%support_folder=[pwd '/support_files'];
addpath(genpath(support_folder));

if run_locally ==1
    %Some hardcodes:
    wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
else
    settings=settings_comparematrices;%
    np=size(settings.path,2);
    warning('off') %supress addpath warnings to nonfolders.
    for i=1:np
        addpath(genpath(settings.path{i}));
    end
    rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
    rmpath('/home/exacloud/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti')
    wb_command=settings.path_wb_c; %path to wb_command
    warning('on')
end

%%Step 1
%Load data
probability_maps_list = importdata(path_to_probability_maps_list);
probability_mat = zeros(91282, size(probability_maps_list,1));
for i=1:size(probability_maps_list,1)
    disp(['Loading ' num2str(i)])
    this_nets_path = probability_maps_list{i};
    cii=ciftiopen(this_nets_path,wb_command);
    this_net_probability = cii.cdata;
    probability_mat(:,i) = this_net_probability;
    
end

%%Step 2 find network of max probability at each grayordinate.

[~, max_vals] = max(probability_mat,[],2);
max_vals = max_vals';
%Step 2.5 adjust network indicies to account for no network 4 or 6.

greater_than_3 = find(max_vals >3);
max_vals(greater_than_3) = max_vals(greater_than_3)+1;

greater_than_5 = find(max_vals >5);
max_vals(greater_than_5) = max_vals(greater_than_5)+1;


%% save dscalar
cii.cdata = max_vals';
if strcmp(output_name(end-11:end),'dscalar.nii') ==0
    ciftisave(cii,[output_folder filesep output_name '.dscalar.nii'],wb_command)
    output_name =[output_name '.dscalar.nii'];
else
    ciftisave(cii,[output_folder filesep output_name],wb_command)
end
if run_clean_up_script ==1
    clean_dscalars_by_size([output_folder filesep output_name],[],[],[],[],30,[],0,1,1)
else
    disp('No clean-up option selected.  Your data may have small islands of networks, if you care about that sort of thing.')
end

disp('Done.')

end