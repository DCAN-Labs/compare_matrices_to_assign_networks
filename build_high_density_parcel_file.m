function parcel = build_high_density_parcel_file(assignments_dscalar_or_vector,output_name)

% R.Hermosillo 12/13/2022
% This function build a parcel file for your individual-specific network
%file.  A parcel, which contains the RGB, and indicies of each ROI in the networks file, can be handy when using the showM function.
%  
% Inputs are are: 
% 1) a dscalar of networks of a vector of assignments.
% 2) an output name.  A string of characters.  '_parcel .mat will be added to the end of the string.'
% 
% Outputs are:
% 1) a .mat file with a struture variable called "parcel" inside.
% 
% NOTE: this code, assumes that the order follows a similar convention to
% the template matching networks. And would need to be modified to use
% different parcel formats.


this_code = which('simple_cifti_average');
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
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/internal/utilities/MergeTimeSeries'));

warning('on')
wb_command=settings.path_wb_c; %path to wb_command

%% start

%gparcel = parcel;
%TM80dense = ciftiopen('/home/rando149/shared/projects/ADHD_MedChal/TMprobabilistic80.networks_pergrayordinate.32k_fs_LR.dscalar.nii',wb_command);
if isnumeric(assignments_dscalar_or_vector)
    TM80dscalar = assignments_dscalar_or_vector;
else
    TM80dense = ciftiopen(assignments_dscalar_or_vector,wb_command);
    TM80dscalar = TM80dense.cdata;
end

nets = unique(nonzeros(TM80dscalar));

%load('/home/rando149/shared/projects/Polyvertexscore/HumanGordon_parcel.mat');
if ismember(18,nets) ==1 % check to see if the scan is included.  If not load the cannonical network names.
    load('/home/rando149/shared/projects/Polyvertexscore/parcel_probability_map_wscan.mat','parcel');
else
    load('/home/rando149/shared/projects/Polyvertexscore/parcel_probability_map.mat','parcel');
end

new_parcel = parcel;
TM80dscalar_reduced = nonzeros(TM80dscalar);
for i=1:size(nets,1)+3
    try
        netrow = find([new_parcel.power_val] ==i); % find the network number that is equal
        if isempty(netrow)
            disp(['No net:  ' num2str(i)]);
        else
            greys = find(TM80dscalar_reduced == i);
            new_parcel(netrow).ix=greys;
            new_parcel(netrow).n=size(greys,1);
        end
    catch
        disp(['No net:  ' num2str(i)]);
    end
end
parcel = new_parcel;
%save('/home/faird/shared/projects/AnitaOHSUVAcollab/code/TMprobabilistic80.networks_pergrayordinate.32k_fs_LR_parcel.mat','parcel')
disp('Saving parcel...')
disp([output_name '_parcel.mat'])
save([output_name '_parcel.mat'],'parcel');

disp('Done saving parcel file.')
