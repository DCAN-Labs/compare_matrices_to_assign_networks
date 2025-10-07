function network_encroachment(reference_dscalar,test_dscalar, encrouching_network_number,outputname,save_encroachment_dscalar,make_verbose,contrain_to_lobe,path_to_lobe_dlabel)


% This code will gnerate a map of the regions that the selected network is "encourcing into."
% R. Hermosillo 09/12/2024

% inputs are:
% the full path to a refeernce dscalar
% THe full path to asubject test dscalar
% the network number of the encrouching network.
% the output name of the dscalar (without the dscalar part).

%Note that the "ismember" function only looks at regions in the subject encouching network.  Regions that have "shrunk" will not be observed without modifying the code."

all_labels = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON','SCAN'};
net_names = {'DMN','Vis','FP','', 'DAN','', 'VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON','', 'SCAN'};

%add_label_percentages =1;

%right_labels = {'AUD', 'CO','DAN','DMN','FP','MTL','PMN','PON','SAL', 'SMD','SML' , 'Tpole', 'VAN','VIS'};
%right_labels = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};
possible_net_nums = [1 2 3  5  7 8 9 10 11 12 13 14 15 16 18];
%donut_label_font_size =10;
run_locally=0;
if run_locally ==1
    %Some hardcodes:
    wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
else
    
    this_code = which('template_matching_RH');
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
    addpath(genpath('/projects/standard/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti')) % remove non-working gifti path included with MSCcodebase
    addpath(genpath('/projects/standard/faird/shared/code/internal/utilities/plotting-tools'));
    addpath(genpath('/projects/standard/faird/shared/code/internal/utilities/Zscore_dconn'));
    warning('on')
    if exist('wb_command','var')==1
    else
        wb_command=settings.path_wb_c; %path to wb_command
    end
end

conc = strsplit(test_dscalar, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    test_dscalar = importdata(test_dscalar);
else
    test_dscalar = {test_dscalar};
end

if contrain_to_lobe ==1
    ft_cii=ft_read_cifti_mod(path_to_lobe_dlabel);
    %ft_cii=ft_read_cifti_mod('/panfs/jay/groups/14/znahas/shared/projects/PCS/network_surface_area/encroachment/frontal_cortex_only/HumanLobes_frontal_lobe_only.32k_fs_LR.dlabel.nii');
    label_array=importdata('/panfs/jay/groups/14/znahas/shared/projects/PCS/network_surface_area/encroachment/frontal_cortex_only/HCP_frontal_lobe_ROIs_only.32k_fs_LR_index_list.txt');
    label_values = label_array(:,1);
    label_indices=ismember(ft_cii.data,label_values);
end


%% Load data
Ccii = ciftiopen(reference_dscalar,wb_command);
C = Ccii.cdata;
if contrain_to_lobe ==1
    label_indices_size_adj=zeros(size(C,1),1);
    label_indices_size_adj(1:size(label_indices),1) = label_indices;
    C_temp=zeros(size(C,1),1);
    
    for i = 1:size(C,1)
        if label_indices_size_adj(i) ==1
            C_temp(i)=C(i); % write this to a temporary variable to check the data;
        end
    end
    C= C_temp;
end

ref_net_indices = find(C == encrouching_network_number);
ref_net_indices_cort = find(C(1:59412) == encrouching_network_number); %hardcode warning
if contrain_to_lobe ==0
    ref_net_indices_subc = find(C(59413:91282) == encrouching_network_number); %hardcode warning
end
isref_net_indices_log = C == encrouching_network_number;

for s=1:size(test_dscalar)
    disp(num2str(s))
    disp(test_dscalar{s})
    Dcii = ciftiopen(test_dscalar{s},wb_command);
    D = Dcii.cdata;
    
    % probably don't need this part:
    %     if contrain_to_lobe ==1
    %       D=D(1:size(ft_cii.data,1),1); %make this smaller
    %         D=D(label_indices);
    %     end
    
    sub_net_indices = find(D == encrouching_network_number);
    sub_net_indices_cort = find(D(1:59412) == encrouching_network_number); %hardcode warning
    if contrain_to_lobe ==0
        sub_net_indices_subc = find(D(59413:91282) == encrouching_network_number); %hardcode warning
    end
    issub_net_indices_log = D == encrouching_network_number;
    
    %diff_scalar_log= issub_net_indices_log ~= isref_net_indices_log; %probably don't use
    
    diff_scalar= ~ismember(sub_net_indices,ref_net_indices);
    diff_scalar_cort= ~ismember(sub_net_indices_cort,ref_net_indices_cort);
    if contrain_to_lobe ==0
        diff_scalar_subc= ~ismember(sub_net_indices_subc,ref_net_indices_subc);
    end
    
    non_sal_indices = sub_net_indices(diff_scalar);
    %non_sal_indices_cort = sub_net_indices(diff_scalar_cort);
    %non_sal_indices_subc = sub_net_indices(diff_scalar_subc);
    if save_encroachment_dscalar ==1
        
        full_ismember_dscalar = zeros(size(D,1),1);
        full_ismember_dscalar(non_sal_indices)=1;
        
        output_dscalar = zeros(size(D,1),1);
        for g=1:size(D,1)
            if full_ismember_dscalar(g) ==1
                output_dscalar(g,1)= C(g); % get the reference network assingment
            end
        end
        
        
        Dcii.cdata=output_dscalar;
        ciftisave(Dcii,[outputname '.dscalar.nii'],wb_command)
    end
    j=1;
    
    for n=1:max(possible_net_nums)
        if n==4 || n==6 || n==17
            
        else
            en_ref_net_idx = sum(output_dscalar ==n);
            ref_net_idx = sum(C ==n);
            encrouched_percent_all(s,j) = en_ref_net_idx/ref_net_idx;
            
            en_ref_net_idx_cort = sum(output_dscalar(1:59412,1) ==n);
            ref_net_idx_cort = sum(C(1:59412,1) ==n);
            encrouched_percent_cort(s,j) = en_ref_net_idx_cort/ref_net_idx_cort;
            
            en_ref_net_idx_subc = sum(output_dscalar(59413:91282,1) ==n);
            ref_net_idx_subc = sum(C(59413:91282,1) ==n);
            encrouched_percent_subc(s,j) = en_ref_net_idx_subc/ref_net_idx_subc;
            
            if make_verbose ==1
                disp(['percent encroachment in all greyordinates for net: ' num2str(n) ' ' char(net_names(j)) ' : ' num2str(encrouched_percent_all(s,j)*100) ])
                disp(['percent encroachment in subcortex for net: ' num2str(n) ' ' char(net_names(j)) ' : ' num2str(encrouched_percent_subc(s,j)*100) ])
            end
            disp(['percent encroachment in cortex for net: ' num2str(n) ' ' char(net_names(j)) ' : ' num2str(encrouched_percent_cort(s,j)*100) ])
        end
        j=j+1;
    end
end

save([outputname '.mat'],'encrouched_percent_cort','encrouched_percent_all','encrouched_percent_subc','reference_dscalar','test_dscalar')
disp('Done.')
end





