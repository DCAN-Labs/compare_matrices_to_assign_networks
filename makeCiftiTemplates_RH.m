function makeCiftiTemplates_RH(dt_or_ptseries_conc_file,motion_file,Zscore_regions,power_motion,remove_outliers)
%%% load consensus, subjects, networks
%consen = ft_read_cifti_mod('/data/cn6/allyd/variants/120_colorassn_minsize400_manualconsensus.dtseries.nii');
%consen = ft_read_cifti_mod('/mnt/max/shared/code/internal/utilities/community_detection/fair/120_colorassn_minsize400_manualconsensus.dtseries.nii');

%some hardcodes:
FD_threshold = 0.2;

if Zscore_regions == 1
    L_size = 29696; %hardcode - number of parcellations
    R_size = 29716;
    S_size = 31870;
else
end

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
        if power_motion == 1
            %use power method
            load(B{i})
            allFD = zeros(1,length(motion_data));
            for j = 1:length(motion_data)
                allFD(j) = motion_data{j}.FD_threshold;
            end
            FDidx = find(round(allFD,3) == round(FD_threshold,3));
            FDvec = motion_data{FDidx}.frame_removal;
            FDvec = abs(FDvec-1);
            if exist('remove_outliers','var') == 0 || remove_outliers == 1
                disp('Removal outliers not specified.  It will be performed by default.')
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                stdev_temp_filename =[file_root '_temp.txt'];
                cmd = [wb_command ' -cifti-stats ' subs{i} ' -reduce STDEV > ' stdev_temp_filename];
                system(cmd);
                clear cmd
                disp('waiting 10 seconds for writing of temp file before reading. Line:258 (not an error)')
                pause(10);  % to give time to write temp file, failed to load with pause(1)
                STDEV_file=load(stdev_temp_filename); % load stdev of .nii file.
                system(['rm -f ' stdev_temp_filename]); %clean up                
                FDvec_keep_idx = find(FDvec==1); %find the kept frames from the FD mask
                Outlier_file=isthisanoutlier(STDEV_file(FDvec_keep_idx),'median'); %find outlier
                Outlier_idx=find(Outlier_file==1); %find outlier indices
                FDvec(FDvec_keep_idx(Outlier_idx))=0; %set outliers to zero within FDvec
                clear STDEV_file FDvec_keep_idx Outlier_file Outlier_idx
                
            else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
            tmask = FDvec;
        else
            tmask = load(B{i});
        end     
        
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

if Zscore_regions == 1
    disp('Converting to Zscores')
    for i=1:length(seedmapsTR)
        for j=1:size(corrs)
            seedmaps{i}(1:29696,j) = zscore(seedmaps{i}(1:29696,j));
            seedmaps{i}(29697:59412,j) = zscore(seedmaps{i}(29697:59412,j));
            seedmaps{i}(59413:91282,j) = zscore(seedmaps{i}(59413:91282,j));
        end
    end
save(['seedmaps_withinregionZscores' file_root_no_ext '.mat'],'seedmapsTR','-v7.3')
else
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
