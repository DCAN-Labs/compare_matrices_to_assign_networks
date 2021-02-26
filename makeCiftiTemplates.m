%%% load consensus, subjects, networks
%consen = ft_read_cifti_mod('/data/cn6/allyd/variants/120_colorassn_minsize400_manualconsensus.dtseries.nii');
consen = ft_read_cifti_mod('/mnt/max/shared/code/internal/utilities/community_detection/fair/120_colorassn_minsize400_manualconsensus.dtseries.nii');
%subs = textread('/data/cn6/allyd/TRsurfaces/allTRlist.txt','%s');
subs = textread('/mnt/max/shared/projects/midnight_scan_club/template_matching/MSC_subjects.txt','%s');
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};
% consen.data=consen.data(1:59412); %%% if surface only


%%% load BOLD data from subjects
%%% extract mean time series for all voxels labeled network j
%%% compute corr between mean and all voxels in that subject
for i=1:length(subs)
    disp(['subject ' subs{i}]);
    
    %%% load subject BOLD data and tmask
    TR = ft_read_cifti_mod(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/' subs{i} '_LR_333_LR_surf_subcort.dtseries.nii']);
    
    tmask = dlmread(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/total_tmask.txt']);
    TR.data = TR.data(:,tmask>0);

    for j=1:length(network_names)
        if j==4 || j==6
            continue
        end
        %disp(['  network ' network_names{j}]);
        inds = consen.data==j;
        subNetAvg= nanmean(TR.data(inds,:),1);
        for voxel=1:length(TR.data);
            goodvox= ~isnan(TR.data(voxel,:));
            corrs(voxel,j)=paircorr_mod(subNetAvg(goodvox)', TR.data(voxel,goodvox)')';
        end  
        clear inds
    end
    clear TR tmask
    seedmapsTR{i} = corrs;
    clear corrs subNetAvg
end


%%% average within networks across subjects
%%% fisher transform for each subject before making average
%%% reverse fisher transform after average
for j=1:length(network_names)
    if j==4 || j==6
        continue
    end
    
    grpNetAve= zeros(size(seedmapsTR{1}),1);
    for i=1:length(subs)
        grpNetAve=grpNetAve + atanh(seedmapsTR{i}(:,j));
    end
    grpNetAve=grpNetAve ./ length(subs);  
    avgSeedmaps{j}=inverseFisherTransform(grpNetAve);
    temp=ft_read_cifti_mod('/data/cn6/allyd/eroded_WB_mask.dtseries.nii');
    temp.data=avgSeedmaps{j};

    surface=temp.data(1:59412); temp.data(consen.data==0)=0; temp.data(1:59412)=surface; %w subcortex
    ft_write_cifti_mod(['/data/cn6/allyd/cifti_TEST_RevisedTemplate_' network_names{j} '_network_surfOnly.dtseries.nii'], temp);
    
    clear grpNetAve temp
end

