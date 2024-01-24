function makeCiftiTemplates_RH(dt_or_ptseries_conc_file,TR,all_motion_conc_file,project_dir,Zscore_regions,power_motion,remove_outliers, surface_only,use_only_subjects_that_pass_motion_criteria,combined_outliermask_provided,include_scan_net)
%%% load consensus, subjects, networks
%consen = ft_read_cifti_mod('/data/cn6/allyd/variants/120_colorassn_minsize400_manualconsensus.dtseries.nii');
%consen = ft_read_cifti_mod('/mnt/max/shared/code/internal/utilities/community_detection/fair/120_colorassn_minsize400_manualconsensus.dtseries.nii');

%some hardcodes:
FD_threshold = 0.2; FD_column = 21;
%FD_threshold = 0.3; %for infant
check_motion_first =1;
minutes_to_use =10;

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
%addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti/fileio/'))
%addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/read_write_cifti/gifti')); % add  new working gifti path included with MSCcodebase
addpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/community_detection/fair/supporting_scripts')
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
%addpath(genpath('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_scripts'))
addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/')); %Add top level folder to get paircorr_mod.m
warning('on')

wb_command=settings.path_wb_c; %path to wb_command
if include_scan_net ==1
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'     ''    'SCAN'};
else
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};
end
%% Check timeseries and motion files
conc = strsplit(dt_or_ptseries_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    subs = importdata(dt_or_ptseries_conc_file);
else
    subs = {dt_or_ptseries_conc_file};
end

subsfound=0;
for i = 1:length(subs)
    if exist(subs{i},'file') == 0
        NOTE = ['Subject Series ' num2str(i) ' does not exist']
        disp(num2str(subs{i}))
        return
    else
        subsfound=subsfound+1;
    end
end
disp([num2str(subsfound) ' of ' num2str(length(subs)) ' timeseries files found.  All series files exist continuing ...'])

if strcmp('conc',conc) == 1
    B = importdata(all_motion_conc_file);
else
    B = {all_motion_conc_file};
end

%Make sure length of motion files matches length of timeseries
if length(subs)==length(B)
else
    disp('length of motion conc file does not match legnth of series con file')
end

% Check that all motion files in conc file exist
subsfound_motion=0;
for i = 1:length(B)
    if exist(B{i},'file') == 0
        NOTE = ['motion file ' num2str(i) ' does not exist']
        disp(num2str(B{i}))
        %return
    else
        subsfound_motion=subsfound_motion+1;
    end
end
disp([num2str(subsfound_motion) ' of ' num2str(length(B)) ' motion files found. All motion files exist continuing ...'])

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
if include_scan_net ==1
%consen = ft_read_cifti_mod('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned_wHCPscan.dscalar.nii');    
consen = ft_read_cifti_mod('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned_wABCDscan.dscalar.nii');    
else
consen = ft_read_cifti_mod(settings.path{6}); %path to dscalar with template labels.
end

%subs = textread('/data/cn6/allyd/TRsurfaces/allTRlist.txt','%s');
%subs = textread('/mnt/max/shared/projects/midnight_scan_club/template_matching/MSC_subjects.txt','%s');
%consen.data=consen.data(1:59412); %%% if surface only
if check_motion_first ==1
    % check_twins_motion(all_motion_conc_file, TR_or_TR_conc, output_dir, output_name, FD_column,FD_conc_file,get_mean_FD,use_outlierdetection_mask_if_possible,split_motion,splits)
[~,output_name] = fileparts(all_motion_conc_file);
   %[all_subjects_minutes,all_mean_FD] = check_twins_motion(all_motion_conc_file, TR, FD_column,[],0,1);
    [all_subjects_minutes,all_mean_FD] = check_twins_motion(all_motion_conc_file, TR, project_dir, output_name ,FD_column,[],0,1,0,0);

         good_subs_idx = find(minutes_to_use<=all_subjects_minutes); %use these subjects
         bad_subs_idx = find(minutes_to_use>all_subjects_minutes); %dont' use these subjects
         really_bad_subs_idx= find(0.5>=all_subjects_minutes); %these subjects have very little data.
end

%Going Forward, only use timeseries and motion from subjects that pass
%minimum number of minutes requirement.
if use_only_subjects_that_pass_motion_criteria ==1
subs = subs(good_subs_idx);
B = B(good_subs_idx);
else
end

%uncomment to plot
% figure()
% subplot(2,1,1)
% histogram(all_subjects_minutes,20);
% xlim([0 max(all_subjects_minutes)+1 ])
% subplot(2,1,2)
% histogram(all_subjects_minutes(good_subs_idx),20);
% xlim([0 max(all_subjects_minutes)+1 ])

%%% load BOLD data from subjects
%%% extract mean time series for all voxels labeled network j
%%% compute corr between mean and all voxels in that subject
if exist([project_dir filesep 'seedmaps_' file_root_no_ext '.mat'],'file') == 2
    disp('loading previously generated data from each subject for template');
    load([project_dir filesep 'seedmaps_' file_root_no_ext '.mat'])
else
    for i=1:length(subs)
        tic
        %for    i=good_subs
        disp(num2str(i));
        disp(['subject ' subs{i}]);
        
        %%% load subject BOLD data and tmask
        %TR = ft_read_cifti_mod(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/' subs{i} '_LR_333_LR_surf_subcort.dtseries.nii']);
        %cii=ciftiopen(subs{i},wb_command);
        %newcii=cii;
        %TR_file = ciftiopen(subs{i}, wb_command);
        %TR = TR_file.cdata;
        timeseries = ft_read_cifti_mod(subs{i});
        %TR=single(cii.data);
        %tmask = dlmread(['/data/cn6/allyd/TRsurfaces/ciftiFiles_TR/' subs{i} '/total_tmask.txt']);
        
        if power_motion == 1
            %use power method
            load(B{i})
            allFD = zeros(1,length(motion_data));
            for motion_colum = 1:length(motion_data)
                allFD(motion_colum) = motion_data{motion_colum}.FD_threshold;
            end
            FDidx = find(round(allFD,3) == round(FD_threshold,3));
            if combined_outliermask_provided ==1
              FDvec = motion_data{FDidx}.combined_removal;  
            else
            FDvec = motion_data{FDidx}.frame_removal;
            end
            FDvec = abs(FDvec-1); % convert "1=remove to 1=keep"
            allmasks_before_outliers_removed_FD02{i} = FDvec;
            
            if combined_outliermask_provided ==0
                if exist('remove_outliers','var') == 0 || remove_outliers == 1
                    
                    disp('Removal outliers not specified.  It will be performed by default.')
                    %% additional frame removal based on Outliers command: isoutlier with "median" method.
                    stdev_temp_filename =[char(file_root(1)) '_temp.txt'];
                    addpath('/mnt/max/shared/code/internal/utilities/CensorBOLDoutliers/')
                    [FDvec]= CensorBOLDoutliers(wb_command, subs, i, stdev_temp_filename, FDvec);
                    
                else exist('remove_outliers','var') == 1 && remove_outliers == 0;
                    disp('Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed. Maybe you have already performed outlier detection.');
                end
                
                good_frames_idx = find(FDvec == 1);
                good_minutes = (length(good_frames_idx)*TR)/60; % number of good minutes in your data
                
                %disp(['Available possible minutes is ' good_minutes'])
                
                if good_minutes < 0.5 % if there is less than 30 seconds, don't generate the correlation matrix
                    subject_has_too_few_frames = ['Subject ' num2str(i) ' has less than 30 seconds of good data']
                    bad_mover(m) = i;
                elseif minutes_to_use > good_minutes % if there is not enough data for your subject, just generate the matrix with all available frames
                    fileID = fopen([char(project_dir) filesep char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_All_Good_Frames.txt'],'w');
                    fprintf(fileID,'%1.0f\n',FDvec);
                    fclose(fileID);
                    
                elseif minutes_to_use <= good_minutes % if there is enough data, match the amount of data used for your subject matrix to the minutes_limit
                    good_frames_needed = round(minutes_to_use*60/TR); %number of good frames to randomly pull
                    rand_good_frames = sort(randperm(length(good_frames_idx),good_frames_needed));
                    FDvec_cut = zeros(length(FDvec),1);
                    ones_idx = good_frames_idx(rand_good_frames);
                    FDvec_cut(ones_idx) = 1; % the new vector that should match good frames with the minutes limit
                    
                    fileID = fopen([char(project_dir) filesep char(orig_motion_filename) '_' num2str(FD_threshold) '_cifti_censor_FD_vector_' num2str(minutes_limit) '_minutes_of_data_at_' num2str(FD_threshold) '_threshold.txt'],'w');
                    fprintf(fileID,'%1.0f\n',FDvec_cut);
                    fclose(fileID);
                    
                end
            else
                disp('Frames with outliers in BOLD std dev not removed. Maybe you have already performed outlier detection.');
            end
            tmask = logical(FDvec);
            %tmask = FDvec;
        else
            tmask = logical(load(B{i})); %load .txt file that has a list of all masks.
            %tmask = load(B{i}); %load .txt file that has a list of all masks.
            
        end
        
        allmasks_outliers_removed_FD02{i} = FDvec; %save for later
       
        timeseries_temp = timeseries;
        %clear timeseries
        %timeseries.data = timeseries_temp(:,tmask>0);
        %timeseries.data = timeseries_temp(:,tmask);
        %timeseries.data = timeseries.data(:,tmask>0); % censor time series with mask.
                timeseries.data = timeseries.data(:,tmask); 
        for j=1:length(network_names)
            if j==4 || j==6 || j==17
                continue
            end
            %disp(['  network ' network_names{j}]);
            inds = consen.data==j;
            subNetAvg= nanmean(timeseries.data(inds,:),1);
            if surface_only ==1
                for voxel=1:length(timeseries.data(1:ncortgrey))
                    goodvox= ~isnan(timeseries.data(voxel,:));
                    corrs(voxel,j)=paircorr_mod(subNetAvg(goodvox)', timeseries.data(voxel,goodvox)')';
                end
            else
                for voxel=1:length(timeseries.data)
                    goodvox= ~isnan(timeseries.data(voxel,:));
                    corrs(voxel,j)=paircorr_mod(subNetAvg(goodvox)', timeseries.data(voxel,goodvox)')';
                end
            end
            clear inds
        end
        clear timeseries tmask
        seedmapstimeseries{i} = corrs;
        clear subNetAvg cii
        toc
    end
    save(['seedmaps_' file_root_no_ext '.mat'],'seedmapstimeseries','allmasks_outliers_removed_FD02','allmasks_before_outliers_removed_FD02', '-v7.3')
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
        if sum(isnan(seedmapstimeseries{i}(:,j))) > 0
            disp(['This subject ' num2str(i) ' has Nans'])
            disp(['File with nans is : ' subs{i} ])
            disp(['number of greyordinates with nans = ' num2str(sum(isnan(seedmapstimeseries{i}(:,j))))])
            badsubidx(k,1) = i; k = k+1;
        else
        end
    end
end

if exist('badsubidx','var') == 0
    disp('Congratulations, your data has no nans.');
    badsubidx = [];
else
    badsubidx = unique(badsubidx);
    cleansubs(badsubidx,:) = [];
    %subs = cleansubs;
    cleanseedmapstimeseries = seedmapstimeseries;
    cleanseedmapstimeseries(badsubidx) = [];
    seedmapstimeseries = cleanseedmapstimeseries;
    disp([num2str(length(badsubidx)) ' subjects were removed from average for "Nans" in greyordinates.'])
end

%% Start averaging
for j=1:length(network_names)
    if j==4 || j==6
        continue
    end
    grpNetAve= zeros(size(seedmapstimeseries{1},1),1);
    for i=1:length(cleansubs)
        grpNetAve=grpNetAve + atanh(seedmapstimeseries{i}(:,j));
    end
    grpNetAve=grpNetAve ./ length(cleansubs);
    %avgSeedmaps{j}=inverseFisherTransform(grpNetAve);
    avgSeedmaps{j}=tanh(grpNetAve);
    
    if surface_only == 0
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
    
    temp_file=ciftiopen(cleansubs{i},wb_command);
    if surface_only ==1
        %temp_file=ciftiopen('/mnt/max/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii',wb_command);
        temp_file= ciftiopen(settings.path{10},wb_command); % path to a surface_only dtseries.nii
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
        %ft_write_cifti_mod([project_dir filesep '/support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_networkZscored.dtseries.nii'], temp);
        ciftisave(temp_file,[project_dir filesep 'seedmaps_' file_root_no_ext '_' network_names{j} '_networkZscored.dtseries.nii'],wb_command);
        seed_matrix(:,j) =avgSeedmaps{j};
    else
        %ft_write_cifti_mod([project_dir filesep 'support_files/seedmaps_' file_root_no_ext '_' network_names{j} '_network.dtseries.nii'], temp);
        ciftisave(temp_file,[project_dir filesep 'seedmaps_' file_root_no_ext '_' network_names{j} '_network.dtseries.nii'],wb_command);
        seed_matrix(:,j) =avgSeedmaps{j};
    end
    
    clear grpNetAve temp
end
%save all maps
if Zscore_regions == 1
    save([project_dir filesep 'seedmaps_' file_root_no_ext '_all_networksZscored.mat'],'B','seed_matrix','subs','Zscore_regions','power_motion','remove_outliers','surface_only','use_only_subjects_that_pass_motion_criteria','combined_outliermask_provided','include_scan_net','bad_subs_idx','good_subs_idx','cleansubs');
else
    save([project_dir filesep 'seedmaps_' file_root_no_ext '_all_networks.mat'],'B','seed_matrix','subs','Zscore_regions','power_motion','remove_outliers','surface_only','use_only_subjects_that_pass_motion_criteria','combined_outliermask_provided','include_scan_net','bad_subs_idx','good_subs_idx','cleansubs');
end

disp('Done making network templates based on the subjects you provided.');
end
