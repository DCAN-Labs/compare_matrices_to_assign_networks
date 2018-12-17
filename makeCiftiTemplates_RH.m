function makeCiftiTemplates_RH(dt_or_ptseries_conc_file,motion_file)
%%% load consensus, subjects, networks
%consen = ft_read_cifti_mod('/data/cn6/allyd/variants/120_colorassn_minsize400_manualconsensus.dtseries.nii');
%consen = ft_read_cifti_mod('/mnt/max/shared/code/internal/utilities/community_detection/fair/120_colorassn_minsize400_manualconsensus.dtseries.nii');

%% Adding paths for this function
support_folder=[pwd '/support_files'];
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end

rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
warning('on')

wb_command=settings.path_wb_c; %path to wb_command
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};

%% Check timeseries and motion files
conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    subs = importdata(dt_or_ptseries_conc_file);
else
    subs = {dt_or_ptseries_conc_file};
end

for i = 1:length(subs)
    if exist(subs{i}) == 0
        NOTE = ['Subject Series ' num2str(i) ' does not exist']
        return
    else
    end
end
disp('All series files exist continuing ...')

if strcmp('conc',conc) == 1
    B = importdata(motion_file);
else
    B = {motion_file};
end

%Make sure length of motion files matches length of timeseries
if length(subs)==length(B)
else
    disp('length of motion conc file does not match legnth of series con file')
end

% Check that all motion files in conc file exist
for i = 1:length(B)
    if exist(B{i}) == 0
        NOTE = ['motion file ' num2str(i) ' does not exist']
        return
    else
    end
end
disp('All motion files exist continuing ...')

file_dir = dt_or_ptseries_conc_file;
file_split = strsplit(file_dir,'/');
file_name = char(file_split(end));
file_root = strsplit(file_name,'.');
file_root_no_ext = char(file_root(1));
 %template_cii=ciftiopen(settings.path{6}, wb_command); %parcellated networks labels path (pscalar).
 %template_labels = template_cii.cdata;

consen = ft_read_cifti_mod(settings.path{6}); %path to dscalar with template labels.
%subs = textread('/data/cn6/allyd/TRsurfaces/allTRlist.txt','%s');
%subs = textread('/mnt/max/shared/projects/midnight_scan_club/template_matching/MSC_subjects.txt','%s');
% consen.data=consen.data(1:59412); %%% if surface only


%%% load BOLD data from subjects
%%% extract mean time series for all voxels labeled network j
%%% compute corr between mean and all voxels in that subject
if exist([pwd '/seedmaps_' file_root_no_ext '.mat'],'file') == 2
    disp('loading previously generated data from each subject for template');
    load([pwd '/seedmaps_' file_root_no_ext '.mat'])
else
    for i=1:length(subs)
        disp(['subject ' subs{i}]);
        
        %%% load subject BOLD data and tmask
        %TR = ft_read_cifti_mod(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/' subs{i} '_LR_333_LR_surf_subcort.dtseries.nii']);
        %cii=ciftiopen(subs{i},wb_command);
        %newcii=cii;
        TR = ft_read_cifti_mod(subs{i});
        %TR=single(cii.data);
        
        %tmask = dlmread(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/total_tmask.txt']);
        %TR.data = TR.data(:,tmask>0);
        tmask = load(B{i});
        
        
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
        clear corrs subNetAvg cii
    end
    save(['seedmaps_' file_root_no_ext '.mat'],'seedmapsTR','-v7.3')
end

%%% average within networks across subjects
%%% fisher transform for each subject before making average
%%% reverse fisher transform after average
if exist([pwd '/support_files/seedmaps_' file_root_no_ext '_all_networks_fromsurfonly.mat'],'file') == 2
    disp('loading previously generated data for template');
    load([pwd '/support_files/seedmaps_' file_root_no_ext '_all_networks_fromsurfonly.mat'])
else
for j=1:length(network_names)
    if j==4 || j==6
        continue
    end
    
    grpNetAve= zeros(size(seedmapsTR{1},1),1);
    for i=1:length(subs)
        grpNetAve=grpNetAve + atanh(seedmapsTR{i}(:,j));
    end
    grpNetAve=grpNetAve ./ length(subs);  
    %avgSeedmaps{j}=inverseFisherTransform(grpNetAve);
    avgSeedmaps{j}=tanh(grpNetAve);
    
    temp=ft_read_cifti_mod(subs{i});
    temp.data=avgSeedmaps{j};

    surface=temp.data(1:59412); temp.data(consen.data==0)=0; temp.data(1:59412)=surface; %w subcortex
    %RJH: not sure what the purpose is of ther pervious line.  Possibly
    %made for a tempalte that is mising data?
    %ft_write_cifti_mod(['/data/cn6/allyd/cifti_TEST_RevisedTemplate_' network_names{j} '_network_surfOnly.dtseries.nii'], temp);
    ft_write_cifti_mod([pwd '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_network_surfOnly.dtseries.nii'], temp);
    seed_matrix(:,j) =temp.data;
    save([pwd '/support_files/seedmaps_' file_root_no_ext '_all_networks_fromsurfonly.mat'],'seed_matrix');
    clear grpNetAve temp
end
end
disp('Done making network templates based on the subjects you provided.');
end
