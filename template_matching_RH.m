function [new_subject_labels, output_cifti_scalar_name] = template_matching_RH(dconn_filename, data_type, template_path,transform_data,output_cifti_name, cifti_output_folder, wb_command, make_cifti_from_results,allow_overlap,overlap_method,surface_only,already_surface_only)

%subjectlist = subject (e.g. dconn.nii)
%data_type = "parcellated" or "dense" connectivity matrix
%template_path = path to .mat file that has th network templates.
%transform_data =  if you want to convert you data before comparing to your template, use can use 1 of 2 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearons' or 'Convert_to_Zscores' or use no tranformation
%output_cifti_name = output_cifti_name pretty clear
%cifti_output_folder = your project directory
%wb_command = path to run HCP workbench command.
%make_cifti_from_results = set to 1 if you want to save your results as a cifti. 0 will not save anything.
%allow_overlap = set to 1 if you're using overlapping networks in your cifti (Your input networks file will likely be a .dtseries.nii)
%overlap_method =  currently, the only supported method is
%"smooth_then_derivative"

%%% load subjects, network info %%%
%load /data/cn6/allyd/kelleyData/pnm_list_revised.mat
%load /data/cn6/allyd/kelleyData/allSubjects.mat %%% allSubjects
%load /data/cn6/allyd/kelleyData/subjectswith20min.mat %%% subs_25min


%allow_overlap = 1;
%overlap_method = 'smooth_then_derivative';
%thresholds = [1:0.25:3.5];
Zscore_eta = 0; % not necessary to zscore
%% Adding paths for this function
this_code = which('template_matching_RH');
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
%rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
%rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/internal/utilities/Zscore_dconn'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command

network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};

if isnumeric(make_cifti_from_results)==1
else
    if strcmp(make_cifti_from_results,'none') == 1
    else
        make_cifti_from_results = str2double(make_cifti_from_results);
    end
end

if isnumeric(allow_overlap)==1
else
        allow_overlap = str2double(allow_overlap);
end

if isnumeric(surface_only)==1
else
        surface_only = str2double(surface_only);
end

if isnumeric(already_surface_only)==1
else
        already_surface_only = str2double(already_surface_only);
end

switch data_type
    case 'parcellated'
        disp('data_type is parcellated')
    case 'dense'
        disp('data_type is dense')
    otherwise
        disp('data_type must be parellated or dense')
        return
end

disp('NOTE: Your template-matching threshold should match the units of your template');
%%% threshold template values @ min value %%%
switch transform_data
    case 'Convert_FisherZ_to_r'
        TEMPLATEMINIMUM = 0.37;
    case 'Convert_r_to_Fisher'
        TEMPLATEMINIMUM = 0.37;
    case 'Convert_to_Zscores'
        TEMPLATEMINIMUM = 1.00;
    otherwise
        TEMPLATEMINIMUM = 0.37;
end

disp(['template minimum is set at ' num2str(TEMPLATEMINIMUM)]);



%load('/data/cn6/allyd/TRsurfaces/templatemat_fullCortex.mat'); %template matrix
load(template_path);
cifti_template_mat_full =seed_matrix;
cifti_template_mat_full(cifti_template_mat_full<= TEMPLATEMINIMUM) = nan;

if iscell(dconn_filename) ==1
else
    dconn_filename = {dconn_filename};
end

%change name for Zscored data
if strcmp(transform_data,'Convert_to_Zscores') ==1
    output_cifti_name = [output_cifti_name '_Zscored'];
else
end

for sub = 1:length(dconn_filename)
    
    if exist([cifti_output_folder '/' output_cifti_name '.dscalar.nii' ], 'file') == 2
        disp('Template_matching dscalar already found for this subject. loading...')
        subject_cii =ciftiopen([cifti_output_folder '/' output_cifti_name '.dscalar.nii'], wb_command);
        new_subject_labels = subject_cii.cdata;
        
        switch data_type
            case 'parcellated'
                output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.pscalar.nii' ];
            case 'dense'
                output_cifti_scalar_name  = [cifti_output_folder '/' output_cifti_name '.dscalar.nii' ];
        end
        
    else
        
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
            tic
            if surface_only ==1
                if already_surface_only == 1
                    %do nothing
                else % take out subcortical connections from dconn.
                    large_subjectdconn = char(dconn_filename);
                    [dconn_filename] = surface_only_dconn(char(dconn_filename),'inferred');
                    dconn_filename = [dconn_filename '.dconn.nii'];
                    disp('Removing connectivity matrix  since a smaller on eiwth surface only has been saved.')
                    %cmd = (['rm -f ' num2str(large_subjectdconn)]);
                    %system(cmd)
                end
            else
            end
            
            
            disp('opening subject dconn...')
            switch transform_data
                case 'Convert_FisherZ_to_r'
                    subject_cii=ciftiopen(char(dconn_filename), wb_command); %dconn path
                    corr_mat_full = single(subject_cii.cdata);
                    
                    if range(corr_mat_full)>2
                        disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
                    elseif range(corr_mat_full) <=2
                        disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
                    end
                    disp('Converting from Fisher Z to Person (tanh).')
                    corr_mat_full = tanh(corr_mat_full);
                    
                case 'Convert_r_to_Fisher'
                    subject_cii=ciftiopen(char(dconn_filename), wb_command); %dconn path
                    corr_mat_full = single(subject_cii.cdata);
                    if range(corr_mat_full)>2
                        disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
                    elseif range(corr_mat_full) <=2
                        disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
                    end
                    disp('Converting from Person to Fisher Z(atanh).')
                    corr_mat_full = atanh(corr_mat_full);
                    
                case 'Convert_to_Zscores'
                    disp('Converting from to Z-scores region-wise. Input dconn can either be pearson or fisherZ.')
                    %                 addpath(genpath('/mnt/max/shared/code/internal/utilities/Zscore_dconn/'))
                    %                 addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/Zscore_dconn'))
                    if surface_only ==1
                        Zdconn = Zscore_dconn_surface_only(char(dconn_filename{sub}),'inferred');
                    else
                        Zdconn = Zscore_dconn(char(dconn_filename{sub}),'inferred');
                    end
                    disp(['loading Zscored dconn: ' char(Zdconn)])
                    subject_cii=ciftiopen([char(Zdconn)], wb_command); %dconn path
                    corr_mat_full = single(subject_cii.cdata);
                    %[~,output_cifti_name] = fileparts(output_cifti_name);
                    %output_cifti_name = [B '.nii'];
                    
                otherwise
                    disp('Data transformation method not found. No tranformation will be applied.  If tranformation is desired, please select: "Convert_FisherZ_to_r" or "Convert_r_to_Pearons" or "Convert_to_Zscores".')
                    subject_cii=ciftiopen(char(dconn_filename), wb_command); %dconn path
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
            if exist('allow_overlap','var') == 1 && allow_overlap == 1
                if allow_overlap == 1
                    disp('Calculating overlap')
                    MuI_threshhold_all_networks = findoverlapthreshold(eta_to_template_vox,network_names,Zscore_eta, overlap_method);
                else
                end
            else
            end
            [x, new_subject_labels] = max(eta_to_template_vox,[],2); %find max for template matching
            
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
            disp(['Saving .mat file: ' cifti_output_folder '/' output_cifti_name '.mat'])
             save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','new_subject_labels','network_names','-v7.3')
%             new_subject_labels = eta_subject_index;
            
            switch transform_data
                case 'Convert_to_Zscores'
                    unix(['rm -f ' char(Zdconn) ])
                otherwise
            end
            
            
            % if make_cifi_from_results == 1; %DF: This variable was not defined, I changed to make a cifti if the mat file exists.
            
            %if exist([cifti_output_folder '/' output_cifti_name '.mat'],'file') == 2
            disp('saving file to cifti')
            if surface_only ==1
                saving_template =ciftiopen(settings.path{12}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
            else %assume 91282
                saving_template =ciftiopen(settings.path{8}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
            end
            
            saving_template.cdata = single(new_subject_labels);
            
            %addpath('/mnt/max/shared/code/internal/utilities/corr_pt_dt/support_files');
            disp('Saving new scalar')
            save(saving_template, [cifti_output_folder '/' output_cifti_name '.gii'],'ExternalFileBinary')
            %save(saving_template, output_cifti_scalar_name, wb_command)
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
            
            if exist('allow_overlap','var') == 1 && allow_overlap == 1
                %open example dtseries.nii
                
                disp('Saving overlap files as dtseries')
                if surface_only ==1
                    cii =ciftiopen(settings.path{11}, wb_command);
                else %assume 91282
                    
                    cii =ciftiopen(settings.path{10}, wb_command);
                end
                
                switch overlap_method
                    case'hist_localmin'
                        for j=1:length(network_names)
                            if j==4 || j ==6
                                continue
                            end
                            
                            %etaZ(:,j) = zscore(eta_to_template_vox(:,j));
                            etaZabovethreshold = eta_to_template_vox(:,j) > MuI_threshhold_all_networks(j);
                            overlapeta(:,j) = etaZabovethreshold*j;
                        end
                        series = overlapeta;
                        cii.cdata = uint8(series);
                        ciftisave(cii,[cifti_output_folder '/' output_cifti_name '_overlap_' overlap_method '.dtseries.nii'],wb_command);
                        
                    case 'smooth_then_derivative'
                        for j=1:length(network_names)
                            if j==4 || j ==6
                                continue
                            end
                            
                            %etaZ(:,j) = zscore(eta_to_template_vox(:,j));
                            etaZabovethreshold = eta_to_template_vox(:,j) > MuI_threshhold_all_networks(j);
                            overlapeta(:,j) = etaZabovethreshold*j;
                        end
                        series = overlapeta;
                        cii.cdata = uint8(series);
                        ciftisave(cii,[cifti_output_folder '/' output_cifti_name '_overlap_' overlap_method '.dtseries.nii'],wb_command);
                        
                    case 'etaZ'
                        for k =1:size(thresholds,2)
                            series = squeeze(overlapetaZ(:,:,k));
                            cii.cdata = uint8(series);
                            ciftisave(cii,[cifti_output_folder '/' output_cifti_name '_overlap_threshold_' num2str(thresholds(k)) '.dtseries.nii'],wb_command);
                        end
                    otherwise
                        disp('Overlap network method not specified.')
                end
                % make sure that"assign_unassinged is set to zero" for
                % overlapping networks. Otherwise the whole brain will be assinged to every network.
                clean_dscalars_by_size([cifti_output_folder '/' output_cifti_name '_overlap_' overlap_method '.dtseries.nii'],[],[],[],[],30,[],0,0,0);
            else
            end
            
        end
        % else
        %     disp(['Something went wrong.' cifti_output_folder '/' output_cifti_name '.mat not found.  Data was not made into a cifti.' ]);
        % end
    end
end

disp(['Cleaning: ' output_cifti_scalar_name]);
[outname] = clean_dscalars_by_size(output_cifti_scalar_name,[],[],[],[],30,[],0,1,1); %do not change the make concensus to 1 here. It will ruin your network assingments.
disp(['Clean file: ' outname '.dscalar.nii'])

% cmd = ['mv ' cifti_output_folder '/' output_cifti_scalar_name];
% unix(cmd)
% clear cmd

clear subs nets tmask i
end
