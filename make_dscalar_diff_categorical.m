function [alloutscalar,outscalar] = make_dscalar_diff_categorical(dscalarApath,dscalarBpath,outname, outputmap_type)

%R. Hermosillo 6/31/2019

% This code makes a difference map from two categrical dscalars.  
%If the values at a given greyordinate are the same, then it labels those as greyordinate as a Nan.

%inputs are:
%dscalarApaths = can be a .conc file or the path to a single dscalar.nii
%dscalarBpaths = can be a .conc file or the path to a single dscalar.nii
%outputname = the name of the output dscalar.
%outputmap_type = options are 'difference' or 'common'.
    %'difference' will outputwhere the network assignment for the two dscalars are different.
    %'common' will output where the two dscalars have networks assignments in common.

%outputs are;
%If a single network is provided for dscalarsA and B, the output will be a dscalar showing where the network differences/commonalities are.
%If conc files are provided, then the resulting output dscalar is the
%probability of a different/common network assignment.


%% Adding paths for this function
this_code = which('make_dscalar_diff_categorical');
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
%rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools/showM'))

warning('on')
wb_command=settings.path_wb_c; %path to wb_command

    % Check to see if there 1 subject or a list of subjects in conc file.
    conc = strsplit(dscalarApath, '.');
    conc = char(conc(end));
    if strcmp('conc', conc)
        pathsAlist = importdata(dscalarApath);
        pathsBlist = importdata(dscalarBpath);
    else
        pathsAlist = {dscalarApath};
        pathsBlist = {dscalarBpath};
    end

    % Validate that all paths in .conc file point to real files
    for i = 1:length(pathsAlist)
        if ~exist(pathsAlist{i}, 'file')
            disp(['Subject scalar ' num2str(i) ' does not exist.'])
            return
        end
    end
     for i = 1:length(pathsBlist)
        if ~exist(pathsBlist{i}, 'file')
            disp(['Subject scalar ' num2str(i) ' does not exist.'])
            return
        end
    end   
    disp('All scalar files exist. Continuing ...')
    
%% Start comparison
for i=1:size(pathsAlist,1)
    
    disp(i)
    Acii = ciftiopen(pathsAlist{i},wb_command);
    Adscalar = single(Acii.cdata);
    
    Bcii = ciftiopen(pathsBlist{i},wb_command);
    Bdscalar = single(Bcii.cdata);
    
    
    outscalar = NaN(size(Adscalar,1),1);
     if i ==1
    alloutscalars = zeros(size(outscalar,1),size(pathsAlist,1));
     end
    %same_net = find(Adscalar == Bdscalar);
    switch outputmap_type
        case 'difference'
            diff_net= find(Adscalar ~= Bdscalar);
            outscalar(diff_net)=Bdscalar(diff_net);
            
        case 'common'
            diff_net= find(Adscalar == Bdscalar);
            outscalar(diff_net)=Bdscalar(diff_net);
        otherwise
            error('outputmap_type variable must be set to "difference" or "common" (as strings).')
    end
    Bcii.cdata = outscalar;
    if strcmp('conc', conc)
        alloutscalars(:,i) = outscalar;
    else
        
        alloutscalar = [];
        ciftisave(Bcii,[outname '_' outputmap_type '.dscalar.nii'],wb_command)
    end
    clear outscalar
    
    
end

  if strcmp('conc', conc)
        isnanarray=~isnan(alloutscalars(:,:)); %nan are shown as 0s.
        count_diff = sum(isnanarray,2);
        outscalar = count_diff/size(pathsAlist,1);
        Bcii.cdata = outscalar;
        ciftisave(Bcii,[outname '_' outputmap_type '.dscalar.nii'],wb_command)
  end
  
end