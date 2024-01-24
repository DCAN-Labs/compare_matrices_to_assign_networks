function [output_name ] = simple_cifti_average(cifti_conc,output_name, method, plot_data, TMprobability_modules,gordon_modules,other_modules,check_for_nans,remove_subjects_with_nans,average_of_diff_of_cifti_pairs,cifti_conc2,matrix_color_limits_vector)

%this code simple loads in already built ciftis (one at a time) rom your conc file and averages
%them.

%% Adding paths for this function
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
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/internal/utilities/Nan_checker'));

warning('on')
wb_command=settings.path_wb_c; %path to wb_command

if iscell(cifti_conc)
    ciftis = cifti_conc; % untested
elseif istable(cifti_conc)
    prompt = 'Data is a table. In which column of the table are the ciftis located? Enter a numeric value. [e.g 1-10] : ';
    table_col = input(prompt);
    table_conc = cifti_conc(table_col,:);
    cifti_conc = table2cell(table_conc);
    ciftis = cifti_conc;
else % assume outside text file
    ciftis =importdata(cifti_conc);
    if average_of_diff_of_cifti_pairs ==1
     ciftis2 =importdata(cifti_conc2);
    end
end

tic
%check to make sure that ciftis exist
found_files_total=0; missing_files_total =0;% make a "found files counter"
for i = 1:length(ciftis)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(ciftis{i}, 'file') == 0
        disp(['Error Subject dscalar ' num2str(i) ' does not exist'])
        disp(ciftis{i});
        missing_files_total = missing_files_total+1;
        missing_files_indx(missing_files_total) = i;
        missing_files{missing_files_total} = ciftis{i};
        %return
    else
        found_files_total = found_files_total+1;
        found_files_indx(found_files_total) = i;
        found_files{found_files_total} = ciftis{i};
    end
end

if found_files_total ==length(ciftis)
    disp('All series files exist continuing ...')
else
    disp('WARNING: Not all files were found.')
    disp(['Expected to find: ' num2str(length(ciftis))])
    disp(['Files found: ' num2str(found_files_total)])
    disp(['Files missing: ' num2str(missing_files_total)])
    prompt = 'Continue with only found files? [Y/N]: ';
    str = input(prompt, 's');
    if strcmp(str,'y')==1 || strcmp(str,'Y')==1 || strcmp(str,'yes')==1 || strcmp(str,'YES')==1 || strcmp(str,'Yes')==1 || isempty(str)
        disp('Using only found files.')
        ciftis = found_files;
    else
        return
    end
end



cifti_type = strsplit(ciftis{1}, '.');
cifti_exten = char(cifti_type(end-1));

if strcmp('dtseries',cifti_exten) == 1
     ismatfile = 0;
elseif strcmp('dconn',cifti_exten) == 1
     ismatfile = 0;
elseif strcmp('dscalar',cifti_exten) == 1
     ismatfile = 0;
elseif strcmp('ptseries',cifti_exten) == 1
     ismatfile = 0;
elseif strcmp('pscalar',cifti_exten) == 1
     ismatfile = 0;
elseif strcmp('pconn',cifti_exten) == 1
        ismatfile = 0;
 elseif strcmp('mat',char(cifti_type(end))) == 1 
     %assume that data is d/pconn.
     cifti_exten = 'pconn';
     ismatfile = 1;
else
    
    disp('filetype not supported by nan checker: ')
end

if ismatfile == 1
    
%     load(ciftis{i},'unsorted_reduced_dconn');
%     current_cifti = unsorted_reduced_dconn;
try
    load(ciftis{i},'diff_matrix');
    current_cifti = diff_matrix;
    
    clear unsorted_reduced_dconn
catch
    disp('couldnt find diffdconn variable')
    try
        disp('looking for matrix data...')
        load(ciftis{i},'sorted_dconn1');
       current_cifti = sorted_dconn1;     
    catch
        error('unclear what vairable name contains matrix...');
    end
end

else
    cifti_cii = ciftiopen(ciftis{1},wb_command);
    current_cifti = cifti_cii.cdata;
end

avg_cifti = zeros(size(current_cifti,1),size(current_cifti,2)); % this should be symmtrical
if size(current_cifti,1) ~= size(current_cifti,2)
    disp('Did you know that you correlation matrix is not symmetrical? Maybe youre not using matrices.')
    %return
end

if strcmp('dtseries',cifti_exten) == 1
    switch method
        
        case 'average_cifti'
            for i = 1:length(ciftis)
                disp(i)
                current_gii_obj = ciftiopen(ciftis{i},wb_command);
                current_cifti = current_gii_obj.cdata;
                avg_cifti = avg_cifti + current_cifti;
            end
            disp('Now to simply divide')
            avg_cifti = avg_cifti/length(ciftis);
            if convert_to_pearson ==1
                avg_cifti = tanh(avg_cifti);
            end
            
        case 'mode_cifti'
            for i = 1:length(ciftis)
                disp(i)
                current_gii_obj = ciftiopen(ciftis{i},wb_command);
                current_cifti = current_gii_obj.cdata;
                all_ciftis(:,:,i) = current_cifti(:,:);
            end
            disp('Getting mode.')
            for n = 1:size(all_ciftis,2)
            all_modes_cifti = mode(squeeze(all_ciftis(:,n,:))'); 
            %avg_cifti(n) = all_modes_cifti(n)'; % call it all_ciftis just to make things easy.
            avg_cifti(:,n) = all_modes_cifti; % call it all_ciftis just to make things easy.
            
            end
        otherwise
            error('Method is not supported.')
            
    end
    
else
    switch method
        
        case 'average_cifti'
            for i = 1:length(ciftis)
                disp(i)
                if ismatfile==1
%                     load(ciftis{i},'unsorted_reduced_dconn');
%                     current_cifti = unsorted_reduced_dconn;
%                    clear unsorted_reduced_dconn;
try
     load(ciftis{i},'diff_matrix');
     current_cifti = diff_matrix; 
catch
     load(ciftis{i},'sorted_dconn1');
     current_cifti = sorted_dconn1;                    
end
                else
                    current_gii_obj = ciftiopen(ciftis{i},wb_command);
                    current_cifti = current_gii_obj.cdata;
                end
                
                if check_for_nans ==1
                    hasnan(i,1) =  cifti_nancheck(ciftis{i});
                else
                    disp('Not checking for nans.')
                    hasnan(i,1) = 0;
                end
                if remove_subjects_with_nans ==1 
                    if hasnan(i,1) ==0
                     if average_of_diff_of_cifti_pairs ==1
                          if ismatfile ==1
                          load(ciftis2{i},'unsorted_reduced_dconn');
                          current_cifti2 = unsorted_reduced_dconn;
                          clear unsorted_reduced_dconn
                          else
                            current_gii_obj2 = ciftiopen(ciftis2{i},wb_command);
                            current_cifti2 = current_gii_obj2.cdata;
                          end
                            current_cifti=current_cifti-current_cifti2;
                            clear current_cifti2;
                     end
                    avg_cifti = avg_cifti + current_cifti;
                    else
                        %don't add subject
                        disp('skipping adding subject with nans.')
                    end
                else
                    if hasnan(i,1) ==1
                        warning('Subject has nans and is being added to average.')
                      if average_of_diff_of_cifti_pairs ==1
                      error('Do not add subjet with nans to average of difference, unless both pairs are removed.')
                      end
                    end
                    avg_cifti = avg_cifti + current_cifti;
                end
            end
            disp('Now to simply divide')
            if remove_subjects_with_nans ==1
                avg_cifti = avg_cifti/(length(ciftis)-sum(hasnan));

            else
                avg_cifti = avg_cifti/length(ciftis);
            end
            if exist('convert_to_pearson','var') ==1
                if convert_to_pearson ==1
                avg_cifti = tanh(avg_cifti);
                end
            end
            
        case 'mode_cifti'
            for i = 1:length(ciftis)
                disp(i)
                current_gii_obj = ciftiopen(ciftis{i},wb_command);
                current_cifti = current_gii_obj.cdata;
                all_ciftis(:,i) = current_cifti(:,1);
            end
            disp('Getting mode.')
            disp('warning nan_check has not been implemented.')
            avg_cifti = mode(all_ciftis'); % call it all_ciftis just to make things easy.
            avg_cifti = avg_cifti';
        otherwise
            error('Method is not supported.')
            clear current_cifti
            
    end
    
end
if strcmp('pconn',cifti_exten) == 1 || strcmp('dconn',cifti_exten) == 1
    if plot_data ==1
        if TMprobability_modules ==1
%             modules = importdata('TMprobability_modules.csv');
%             nets_assigns = modules.data;
%             networks = unique(nets_assigns);
%             [sorted_networks, netsortindx ] = sort(nets_assigns);
%             
%             imagesc(avg_cifti(netsortindx,netsortindx))
load('/panfs/jay/groups/6/faird/shared/projects/AnitaOHSUVAcollab/code/TMprobabilistic80.networks_pergrayordinate.32k_fs_LR_parcel.mat','parcel');

        end
        if gordon_modules ==1
            load('/home/rando149/shared/projects/Polyvertexscore/HumanGordon_parcel.mat','parcel')
        end
        
        if isnumeric(other_modules) ~=1
            load(other_modules,'parcel')
        end
        if average_of_diff_of_cifti_pairs ==1
            showM(avg_cifti,'parcel',parcel,'clims',matrix_color_limits_vector,'line_color',[0 0 0],'line_width',0.5,'fs_axis',8,'fig_wide',7,'one_side_labels',1,'fig_tall',8);
            load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Positive-Negative_ColorMap.mat','pos_neg_cmap');
            colormap(pos_neg_cmap);
        else
            showM(avg_cifti,'parcel',parcel,'clims',matrix_color_limits_vector,'line_color',[0 0 0],'line_width',0.5,'fs_axis',8,'fig_wide',7,'one_side_labels',1,'fig_tall',8);
            load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Positive-Negative_ColorMap.mat','pos_neg_cmap');
            colormap(pos_neg_cmap);           
            
            %colormap jet
        end
        disp('saving image...')
        print([output_name method '.png'], '-dpng', '-r300')
        
    end
end

%current_cponn.cdata = avg_cifti;
%ciftisave(current_cponn,avg_cifti_output_name,wb_command);
if ismatfile ==1
    save([output_name '.mat'],'avg_cifti','-v7.3')
else
current_gii_obj.cdata = avg_cifti;
ciftisave(current_gii_obj,output_name,wb_command)
end

