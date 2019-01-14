function [eta_to_template_vox, eta_subject_index] = template_matching_RH(subjectlist, data_type, template_path,transform_data,output_cifti_name, cifti_output_folder ,wb_command)

%subjectlist = subject (e.g. dconn.nii)
%data_type = "parcellated" or "dense" connectivity matrix
%template_path = path to .mat file that has th network templates.
%transform_data =  if you want to convert you data before comparing to your template, use can use 1 of 2 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearons' or 'Convert_to_Zscores' or use no tranformation
%output_cifti_name = output_cifti_name pretty clear
%cifti_output_folder = your project directory
%wb_command = path to run HCP workbench command.


%%% load subjects, network info %%%
%load /data/cn6/allyd/kelleyData/pnm_list_revised.mat
%load /data/cn6/allyd/kelleyData/allSubjects.mat %%% allSubjects
%load /data/cn6/allyd/kelleyData/subjectswith20min.mat %%% subs_25min 
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};


support_folder=['/mnt/max/shared/code/internal/analyses/compare_matrices/support_files'];
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

warning('off') %supress addpath warnings to nonfolders.
disp('Attempting to add neccesaary paths and functions.')
for i=2:np % start at 2 to skip MSCcodebase
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
wb_command=settings.path_wb_c; %path to wb_command
warning('on')




%%% threshold template values @ min value %%%
TEMPLATEMINIMUM = 0.37;
%load('/data/cn6/allyd/TRsurfaces/templatemat_fullCortex.mat'); %template matrix
load(template_path);
cifti_template_mat_full =seed_matrix;
cifti_template_mat_full(cifti_template_mat_full<= TEMPLATEMINIMUM) = nan;

if iscell(subjectlist) ==1
else
    subjectlist = {subjectlist};
end


for i = 1:length(subjectlist)
        if exist([cifti_output_folder '/' output_cifti_name '.mat']) == 2
        disp('.mat file already reated.  loading...');
        load([cifti_output_folder '/' output_cifti_name '.mat']);
    new_subject_labels = eta_subject_index;
        else
    %subnum = subjectlist(i);
    %index = find(allSubjects==subnum); %%% index of good subject in full subjects list %%%
    
    %%% concatenate all scan sessions dtseries for each subject %%%
    %catData=[];  
    %catTmask=[];
%     for session=1:size(probabilistic_network_map_list{index,3},1)
%         sessionID = char(probabilistic_network_map_list{index,3}(session));
%         tempseries = ft_read_cifti_mod(['/data/cn5/selfRegulation/CIFTIs/cifti_timeseries_normalwall/' sessionID '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
%         catData = [catData tempseries.data];
%         tempTmask = dlmread(['/data/cn5/selfRegulation/finalOE/' sessionID '/total_tmask.txt'])';
%         catTmask = [catTmask tempTmask];
%         clear tempseries tempTmask
%     end
%  
%     catData = subjectlist{i};
%     catTmask = motion_masks{i};
%     %%% apply tmask to cifti timeseries %%%
%     catData = catData(:,catTmask>0);
%     catData = catData(1:59412,:); %hardcode warning
%     %clear catTmask session sessionID
%     clear catTmask
%     
%     
     tic
%    
%     %%% make correlation matrix %%%
%     disp(['Processing subject ' num2str(subnum) '...']);
%     disp('     computing correlation matrix')
%     corr_mat_full = paircorr_mod(catData');
%     clear catTemp catData
    
    %%% if excluding nearby voxels %%%
    %     corr_mat_thr = corr_mat_full;
    %     clear corr_mat_full
    %     corr_mat_thr(dmat<20) = nan;
    %     clear dmat

    

    %%% if template-matching using correlation %%%    
    %     for i=1:size(corr_mat_full,1)
    %        goodvox = ~isnan(corr_mat_full(i,:));
    %        sim_to_template_vox(i,:) = paircorr_mod(corr_mat_full(i,goodvox)',cifti_template(:,goodvox)');
    %     end
    

    %cii=ciftiopen('/mnt/max/shared/projects/hcp_community_detection/Evan_test/Cifti_Community_Detection/distmat_creation/EUGEODistancematrix_XYZ_255interhem_unit8.pconn.nii',path_wb_c);
    %template_cii=ciftiopen('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii', wb_command);

    disp('opening subject dconn...')
    switch transform_data
        case 'Convert_FisherZ_to_r'
            subject_cii=ciftiopen(char(subjectlist), wb_command); %dconn path
            corr_mat_full = single(subject_cii.cdata);
            
            if range(corr_mat_full)>2
                disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
            elseif range(corr_mat_full) <=2
                disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
            end
            disp('Converting from Fisher Z to Person (tanh).')
            corr_mat_full = tanh(corr_mat_full);
            
        case 'Convert_r_to_Pearons'
            subject_cii=ciftiopen(char(subjectlist), wb_command); %dconn path
            corr_mat_full = single(subject_cii.cdata);
            if range(corr_mat_full)>2
                disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
            elseif range(corr_mat_full) <=2
                disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
            end
            disp('Converting from Person to Fisher Z(atanh).')
            corr_mat_full = atanh(corr_mat_full);
            
        case'Convert_to_Zscores'
            disp('Converting from to Z-scores region-wise.')
            addpath('/mnt/max/shared/code/internal/utilities/Zscore_dconn/')
            output_cifti_name = Zscore_dconn(char(subjectlist{i}),'inferred');
            subject_cii=ciftiopen(char(output_cifti_name), wb_command); %dconn path
            corr_mat_full = single(subject_cii.cdata);
            
        otherwise      
        disp('Data transformation method not found. No tranformation will be applied.  If tranformation is desired, please select: "Convert_FisherZ_to_r" or "Convert_r_to_Pearons" or "Convert_to_Zscores".')
            subject_cii=ciftiopen(char(subjectlist), wb_command); %dconn path
            corr_mat_full = single(subject_cii.cdata);
            if range(corr_mat_full)>2
                disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
            elseif range(corr_mat_full) <=2
                disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
            end
        
    end
    
    
    clear subject_cii %save memory
            
   disp('Calculating similarity (eta) to template') 
    %%% compute eta similarity value b/w each vertex and template %%%
     eta_to_template_vox = single(zeros(size(corr_mat_full,1),length(network_names)));
    for i=1:size(corr_mat_full,1)
        if rem(i,5000)==0
            disp([' Calculating voxel ' num2str(i)]);toc;
        end
        for j=1:length(network_names)
            if j==4 || j ==6
                continue
            end            
            %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
            %goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(j,:)));
            goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(:,j))');
            cmap = corr_mat_full(i,goodvox)';
            %tmap = cifti_template_mat_full(j,goodvox)';
            tmap = cifti_template_mat_full(goodvox,j);
            Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
            Mwithin = (tmap+cmap)/2;
            SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
            SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
            eta_to_template_vox(i,j) = 1 - SSwithin/SStot;

            clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
        end
    end   
    
    clear corr_mat_full goodvox i temp

    %%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
    [x, eta_subject_index] = max(eta_to_template_vox,[],2);

    %%% if requiring a minimum eta value for assignment %%%
    %     for j=1:size(eta_subject_index,1)
    %         if x(j)<0.15 | isnan(x(j))
    %             eta_subject_index(j)=0;
    %         end
    %     end
    
    
    %save(['/data/cn6/allyd/kelleyData/eta_to_template_full_cortex/' num2str(subnum) '_eta_to_template' num2str(TEMPLATEMINIMUM) '.mat'],'eta_to_template_vox','-v7.3');
    
    
    %clear sim_to_template_vox j eta_to_template_vox x
    
    %template = ft_read_cifti_mod('/data/cn6/allyd/variants/120_consensus_surfOnly.dtseries.nii');
    %template.data = eta_subject_index;

    
    %%% write out a cifti file of the subject's final winner-take-all map %%%
    %ft_write_cifti_mod(['/data/cn6/allyd/kelleyData/wta_maps_full_cortex/eta_to_template_sub' num2str(subnum) '_' num2str(TEMPLATEMINIMUM) 'templates_nomincorr_wta_map.dtseries.nii'], template)
    
    %clear new_full_mat net_subject_ind_full_matrix network_subject_index y template eta_subject_index temp
    
    toc
            save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','eta_subject_index','-v7.3')
            new_subject_labels = eta_subject_index;
        end
end

% if make_cifi_from_results == 1; %DF: This variable was not defined, I changed to make a cifti if the mat file exists.

if exist([cifti_output_folder '/' output_cifti_name '.mat'],'file') == 2
    disp('saving file to cifti')
    saving_template =ciftiopen(settings.path{8}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
    saving_template.cdata = single(new_subject_labels);
    %addpath('/mnt/max/shared/code/internal/utilities/corr_pt_dt/support_files');
    disp('Saving new scalar')
    save(saving_template, [cifti_output_folder '/' output_cifti_name '.gii'],'ExternalFileBinary') %DF: This save not pointing to the right place fix!
    %save(saving_template, output_cifti_scalar_name, wb_command) %RH sifti save fix.
    switch data_type
        case 'parcellated'
            output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.pscalar.nii' ];
        case 'dense'
            output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.dscalar.nii' ];
    end
    disp('Converting scalar .gii to .nii')
    unix([wb_command ' -cifti-convert -from-gifti-ext ' cifti_output_folder '/' output_cifti_name '.gii ' output_cifti_scalar_name ]);
    disp('Removing .gii')
    unix(['rm -f ' cifti_output_folder '/' output_cifti_name '.gii']);
    unix(['rm -f ' cifti_output_folder '/' output_cifti_name '.dat']);
else
    disp(['Something went wrong.' cifti_output_folder '/' output_cifti_name '.mat not found.  Data was not made into a cifti.' ]);
end
% cmd = ['mv ' cifti_output_folder '/' output_cifti_scalar_name];
% unix(cmd)
% clear cmd

clear subs nets tmask i
end
