function normalize_myelinmaps(input_myelin_scalars_conc)

%% Add necessary paths
addpath ('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices')
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'))

this_code = which('normalize_myelinmaps');
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

%% Load concatenated paths (i.e. paths to ciftis)
%check to see if there one subject or a a list of subjects in conc file.
conc = strsplit(input_myelin_scalars_conc, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    A = importdata(input_myelin_scalars_conc);
else
    A = {input_myelin_scalars_conc};
end

tic
%check to make sure that surface files exist
if exist('A','var') == 1
    for i = 1:length(A)
        if rem(i,100)==0
            disp([' Validating file existence ' num2str(i)]);toc;
        end
        if exist(A{i},'file') == 0
            disp(['NOTE = Subject dscalar ' num2str(i) ' does not exist'])
            disp(A{i});
            return
        else
        end
    end
    disp('All template matching dscalars exist continuing ...')
end



for i=1:size(A)
    
cii = ciftiopen(A{i},wb_command);
mmap = cii.cdata;
zmmap = zscore(mmap);
cii.cdata = zmmap;

disp('saving data')
[cii_dir, Q] = fileparts(A{i});
[~, cii_basename] = fileparts(Q);
output_cii_name = [cii_dir filesep cii_basename '_Zscored.dscalar.nii'];
ciftisave(cii,output_cii_name,wb_command)

end
disp('Done Zscoring myelinmaps')
end