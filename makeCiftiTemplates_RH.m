function makeCiftiTemplates_RH(dt_or_ptseries_conc_file,motion_file,Zscore_regions,power_motion,remove_outliers, surface_only)
%%% load consensus, subjects, networks
%consen = ft_read_cifti_mod('/data/cn6/allyd/variants/120_colorassn_minsize400_manualconsensus.dtseries.nii');
%consen = ft_read_cifti_mod('/mnt/max/shared/code/internal/utilities/community_detection/fair/120_colorassn_minsize400_manualconsensus.dtseries.nii');

%some hardcodes:
FD_threshold = 0.2;
%FD_threshold = 0.3; %for infant

if surface_only ==1
    Zscore_regions = 0;
end

if Zscore_regions == 1
    L_size = 29696; %hardcode - number of parcellations
    R_size = 29716;
    S_size = 31870;
else
end
ncortgrey = 59412;

%% Adding paths for this function
this_code = which('makeCiftiTemplates_RH');
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
addpath(genpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti/fileio/'))
%addpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti/fileio')
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti/gifti') % remove non-working gifti path included with MSCcodebase
%addpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti')
addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_scripts')
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti/fileio/');
addpath(genpath('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_scripts'))
addpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/') %Add top level folder to get paircorr_mod.m
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

if surface_only == 1
    file_root_no_ext =  [file_root_no_ext '_SurfOnly'];
end

%template_cii=ciftiopen(settings.path{6}, wb_command); %parcellated networks labels path (pscalar).
%template_labels = template_cii.cdata;

%consen_file = ciftiopen(settings.path{6}, wb_command); %path to dscalar with template labels.
%consen = consen_file.cdata;
consen = ft_read_cifti_mod(settings.path{6}); %path to dscalar with template labels.
%subs = textread('/data/cn6/allyd/TRsurfaces/allTRlist.txt','%s');
%subs = textread('/mnt/max/shared/projects/midnight_scan_club/template_matching/MSC_subjects.txt','%s');
%consen.data=consen.data(1:59412); %%% if surface only


%%% load BOLD data from subjects
%%% extract mean time series for all voxels labeled network j
%%% compute corr between mean and all voxels in that subject
if exist([code_dir '/seedmaps_' file_root_no_ext '.mat'],'file') == 2
    disp('loading previously generated data from each subject for template');
    load([code_dir '/seedmaps_' file_root_no_ext '.mat'])
else
    for i=1:length(subs)
        disp(['subject ' subs{i}]);
        
        %%% load subject BOLD data and tmask
        %TR = ft_read_cifti_mod(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/' subs{i} '_LR_333_LR_surf_subcort.dtseries.nii']);
        %cii=ciftiopen(subs{i},wb_command);
        %newcii=cii;
        %TR_file = ciftiopen(subs{i}, wb_command);
        %TR = TR_file.cdata;
        TR = ft_read_cifti_mod(subs{i});
        %TR=single(cii.data);
        %tmask = dlmread(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/total_tmask.txt']);
        
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
            allmasks_before_outliers_removed_FD02{i} = FDvec;
            if exist('remove_outliers','var') == 0 || remove_outliers == 1
                
                disp('Removal outliers not specified.  It will be performed by default.')
                %% additional frame removal based on Outliers command: isoutlier with "median" method.
                stdev_temp_filename =[char(file_root(1)) '_temp.txt'];
                addpath('/mnt/max/shared/code/internal/utilities/CensorBOLDoutliers/')
                [FDvec]= CensorBOLDoutliers(wb_command, subs, i, stdev_temp_filename, FDvec);
                
            else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed');
            end
            tmask = FDvec;
        else
            tmask = load(B{i}); %load .txt file that has a list of all masks.
        end
        
        allmasks_outliers_removed_FD02{i} = FDvec; %save for later
       
        %timeseries_temp = TR;
        %clear TR
        %TR.data = timeseries_temp(:,tmask>0);
        
        %TR.data = TR.data(:,tmask>0); % censor time series with mask.
        
        for j=1:length(network_names)
            if j==4 || j==6
                continue
            end
            %disp(['  network ' network_names{j}]);
            inds = consen.data==j;
            subNetAvg= nanmean(TR.data(inds,:),1);
            if surface_only ==1
                for voxel=1:length(TR.data(1:ncortgrey))
                    goodvox= ~isnan(TR.data(voxel,:));
                    corrs(voxel,j)=paircorr_mod(subNetAvg(goodvox)', TR.data(voxel,goodvox)')';
                end
            else
                for voxel=1:length(TR.data)
                    goodvox= ~isnan(TR.data(voxel,:));
                    corrs(voxel,j)=paircorr_mod(subNetAvg(goodvox)', TR.data(voxel,goodvox)')';
                end
            end
            clear inds
        end
        clear TR tmask
        seedmapsTR{i} = corrs;
        clear subNetAvg cii
    end
    save(['seedmaps_' file_root_no_ext '.mat'],'seedmapsTR','allmasks_outliers_removed_FD02','allmasks_before_outliers_removed_FD02', '-v7.3')
end

% if Zscore_regions == 1
%     disp('Converting to Zscores')
%     for i=1:length(seedmapsTR)
%         if exist('corrs','var') == 1% if data is previously loaded corrs might be a variable in workspace
%         else
%             corrs = ones(91282,16); %build a dummy variable that is a 91282 by 16 matrix.
%         end
%         try
%             for j=1:size(corrs,2)
%                 seedmapsTR{i}(1:29696,j) = zscore(seedmapsTR{i}(1:29696,j));
%                 seedmapsTR{i}(29697:59412,j) = zscore(seedmapsTR{i}(29697:59412,j));
%                 seedmapsTR{i}(59413:91282,j) = zscore(seedmapsTR{i}(59413:91282,j));
%             end
%         catch
%             disp(num2str(i))
%             return
%         end
%     end
%     save(['seedmaps_withinregionZscores' file_root_no_ext '.mat'],'seedmapsTR','-v7.3')
% else
% end
%% There are two methods for averaging, depending on whether or nothe data is Zscored.

%1) if your data is not Zscored, i.e. you want seedmaps based on
%correlation, the process is as follows:

%%% 1.1: Average within networks across subjects
%%% 1.2: Fisher transform for each subject before making average
%%% 1.3  Reverse fisher transform after average

%2) If your data is Zscored, we simply average those Zscores together.

%     if exist([code_dir '/support_files/seedmaps_' file_root_no_ext '_all_networks.mat'],'file') == 2
%         disp('loading previously generated data for template');
%         load([code_dir '/support_files/seedmaps_' file_root_no_ext '_all_networks.mat'])
%     else
%     end

%% Check for nans
disp('Checking for nans in seeds...')
k=1; %counter for bad subjects
cleansubs = subs;
for i=1:length(subs)
    for j=1:length(network_names)
        if j==4 || j==6
            continue
        end
        %seedmapsTR{i}(:,j)
        if sum(isnan(seedmapsTR{i}(:,j))) > 0
            disp(['This subject ' num2str(i) ' has Nans'])
            disp(['File with nans is : ' subs{i} ])
            disp(['number of greyordinates with nans = ' num2str(sum(isnan(seedmapsTR{i}(:,j))))])
            badsubidx(k,1) = i; k = k+1;
        else
        end
    end
end

if exist('badsubidx','var') == 0
    disp('Congratulations, your data has no nans.')
else
    badsubidx = unique(badsubidx);
    cleansubs(badsubidx,:) = [];
    subs = cleansubs;
    cleanseedmapsTR = seedmapsTR;
    cleanseedmapsTR(badsubidx) = [];
    seedmapsTR = cleanseedmapsTR;
    disp([num2str(length(badsubidx)) ' subjects were removed from average for "Nans" in greyordinates.'])
end

%% Start averaging
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
    
    if surface_only == 1
        if Zscore_regions == 1
            disp('Converting to Zscores')
            %for i=1:length(seedmapsTR)
            %         if exist('corrs','var') == 1% if data is previously loaded corrs might be a variable in workspace
            %         else
            %             corrs = ones(91282,16); %build a dummy variable that is a 91282 by 16 matrix.
            %         end
            try
                %for j=1:size(corrs,2)
                avgSeedmaps{j}(1:29696) = zscore(avgSeedmaps{j}(1:29696));
                avgSeedmaps{j}(29697:59412) = zscore(avgSeedmaps{j}(29697:59412));
                if surface_only == 0
                    
                    avgSeedmaps{j}(59413:91282) = zscore(avgSeedmaps{j}(59413:91282));
                end
                %end
            catch
                disp(num2str(i))
                return
            end
            %end
            %save(['seedmaps_withinregionZscores' file_root_no_ext '.mat'],'avgSeedmaps','-v7.3')
        else
        end
    else
    end
    
    
    
    %temp=ft_read_cifti_mod(subs{i});
    %temp.data=avgSeedmaps{j};
    
    temp_file=ciftiopen(subs{i},wb_command);
    if surface_only ==1
        %temp_file=ciftiopen('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii',wb_command);
        temp_file= ciftiopen(settings.path{10},wb_command);
        reset_temp_file = zeros(size(temp_file.cdata,1),1);
        temp_file.cdata = reset_temp_file;
    else
    end
    temp_file.cdata= avgSeedmaps{j};
    
    
    
    %surface=temp.data(1:59412); 
    %temp.data(consen.data==0)=0; 
    %temp_file.data(1:59412)=surface; %w subcortex
    %RJH: not sure what the purpose is of ther pervious line.  Possibly
    %made for a tempalte that is mising data?
    %ft_write_cifti_mod(['/data/cn6/allyd/cifti_TEST_RevisedTemplate_' network_names{j} '_network_surfOnly.dtseries.nii'], temp);
    if Zscore_regions == 1
        %ft_write_cifti_mod([code_dir '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_networkZscored.dtseries.nii'], temp);
        ciftisave(temp_file,[code_dir '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_networkZscored.dtseries.nii'],wb_command);
        seed_matrix(:,j) =avgSeedmaps{j};
    else
        %ft_write_cifti_mod([code_dir '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_network.dtseries.nii'], temp);
        ciftisave(temp_file,[code_dir '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_network.dtseries.nii'],wb_command);
        seed_matrix(:,j) =avgSeedmaps{j};
    end
    
    clear grpNetAve temp
end
%save all maps
if Zscore_regions == 1
    save([code_dir '/support_files/seedmaps_' file_root_no_ext '_all_networksZscored.mat'],'seed_matrix');
else
    save([code_dir '/support_files/seedmaps_' file_root_no_ext '_all_networks.mat'],'seed_matrix');
end

disp('Done making network templates based on the subjects you provided.');
end
