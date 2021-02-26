function make_dscalar_diff_categorical(rest_dscalar,task_dscalar,outname, outputmap_type)

% This code makes a difference map from two categrical dscalars.  It
% assumes that the first one is rest.  If the values at a given
% greyordinate are are the same, then it labels those as greyordinate as a Nan.

%variables are: 

%% Adding paths for this function
this_code = which('make_dscalar_diff_categorical');
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
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command


%% Start comparison
    Acii = ciftiopen(rest_dscalar,wb_command);
    Adscalar = single(Acii.cdata);
    
    Bcii = ciftiopen(task_dscalar,wb_command);
    Bdscalar = single(Bcii.cdata);
    
    
    outscalar = NaN(size(Adscalar,1),1);
    %same_net = find(Adscalar == Bdscalar);
    switch outputmap_type
        case 'difference'
    diff_net= find(Adscalar ~= Bdscalar);
    outscalar(diff_net)=Bdscalar(diff_net);
 
        case 'common'
    diff_net= find(Adscalar == Bdscalar);
     outscalar(diff_net)=Bdscalar(diff_net);
    end
    Bcii.cdata = outscalar;
    ciftisave(Bcii,[outname '.dscalar.nii'],wb_command)
    
end