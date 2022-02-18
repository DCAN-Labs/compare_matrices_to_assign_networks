%%% load subjects, network info %%%
load /data/cn6/allyd/kelleyData/pnm_list_revised.mat
load /data/cn6/allyd/kelleyData/allSubjects.mat %%% allSubjects
load /data/cn6/allyd/kelleyData/subjectswith20min.mat %%% subs_25min 
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};

%%% threshold template values @ min value %%%
TEMPLATEMINIMUM = 0.37;
load('/data/cn6/allyd/TRsurfaces/templatemat_fullCortex.mat'); %template matrix
cifti_template_mat_full(cifti_template_mat_full<= TEMPLATEMINIMUM) = nan;

for i = 1:length(subs_20min)
    subnum = subs_20min(i);
    index = find(allSubjects==subnum); %%% index of good subject in full subjects list %%%
    
    %%% concatenate all scan sessions dtseries for each subject %%%
    catData=[];  
    catTmask=[];
    for session=1:size(probabilistic_network_map_list{index,3},1)
        sessionID = char(probabilistic_network_map_list{index,3}(session));
        tempseries = ft_read_cifti_mod(['/data/cn5/selfRegulation/CIFTIs/cifti_timeseries_normalwall/' sessionID '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
        catData = [catData tempseries.data];
        tempTmask = dlmread(['/data/cn5/selfRegulation/finalOE/' sessionID '/total_tmask.txt'])';
        catTmask = [catTmask tempTmask];
        clear tempseries tempTmask
    end
  
    %%% apply tmask to cifti timeseries %%%
    catData = catData(:,catTmask>0);
    catData = catData(1:59412,:);
    clear catTmask session sessionID
    
    tic
    
    %%% make correlation matrix %%%
    disp(['Processing subject ' num2str(subnum) '...']);
    disp('     computing correlation matrix')
    corr_mat_full = paircorr_mod(catData');
    clear catTemp catData
    
    %%% if excluding nearby voxels %%%
    %     corr_mat_thr = corr_mat_full;
    %     clear corr_mat_full
    %     corr_mat_thr(dmat<20) = nan;
    %     clear dmat

    disp('     calculating similarity (eta) to template')

    %%% if template-matching using correlation %%%    
    %     for i=1:size(corr_mat_full,1)
    %        goodvox = ~isnan(corr_mat_full(i,:));
    %        sim_to_template_vox(i,:) = paircorr_mod(corr_mat_full(i,goodvox)',cifti_template(:,goodvox)');
    %     end
    
    
    %%% compute eta similarity value b/w each vertex and template %%%
    for i=1:size(corr_mat_full,1)
        if rem(i,5000)==0
            disp(['     Calculating voxel ' num2str(i)]);
        end
        for j=1:length(network_names)
            if j==4 || j ==6
                continue
            end            
            %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
            goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(:)));
            cmap = corr_mat_full(i,goodvox)';
            tmap = cifti_template_mat_full(j,goodvox)';
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
    
    
    save(['/data/cn6/allyd/kelleyData/eta_to_template_full_cortex/' num2str(subnum) '_eta_to_template' num2str(TEMPLATEMINIMUM) '.mat'],'eta_to_template_vox','-v7.3');
    
    clear sim_to_template_vox j eta_to_template_vox x
    
    template = ft_read_cifti_mod('/data/cn6/allyd/variants/120_consensus_surfOnly.dtseries.nii');
    template.data = eta_subject_index;

    
    %%% write out a cifti file of the subject's final winner-take-all map %%%
    ft_write_cifti_mod(['/data/cn6/allyd/kelleyData/wta_maps_full_cortex/eta_to_template_sub' num2str(subnum) '_' num2str(TEMPLATEMINIMUM) 'templates_nomincorr_wta_map.dtseries.nii'], template)
    
    clear new_full_mat net_subject_ind_full_matrix network_subject_index y template eta_subject_index temp
    
    toc
end
  
clear subs nets tmask i
