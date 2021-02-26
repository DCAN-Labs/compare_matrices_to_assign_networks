function [all_overlaps_vec, all_overlaps_vec_surf_only] = calculate_network_overlap_area(overlap_dtseries_conc_file, motion_censor_conc_file,TR_conc,L_surface_file,R_surface_file)

%parameter
subsetmotion=1;
%TR = 2.5;
minutes_cutoff =[1,4,5, 10];
colormap_limits = [0 30000]; % sets the upper limit of the colormap.
check_file_existence =0;

%net_names = '';

plotcols = size(minutes_cutoff,2)+1;

this_code = which('calculate_network_overlap_area');
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
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')
wb_command=settings.path_wb_c; %path to wb_command

conc = strsplit(overlap_dtseries_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    A = importdata(overlap_dtseries_conc_file);
else
    A = {overlap_dtseries_conc_file};
end
if strcmp('conc',conc) == 1
    B = importdata(motion_censor_conc_file);
else
    B = {motion_censor_conc_file};
end

if strcmp('conc',conc) == 1
    TRs = importdata(TR_conc);
else
    TRs = {TR_conc};
end

if strcmp('conc',conc) == 1
    C = importdata(L_surface_file);
else
    C = {L_surface_file};
end

if strcmp('conc',conc) == 1
    D = importdata(R_surface_file);
else
    D = {R_surface_file};
end

if check_file_existence ==1
    tic
    for i = 1:length(A)
        if rem(i,100)==0
            disp([' Validating file existence ' num2str(i)]);toc;
        end
        if exist(A{i},'file') == 0
            disp(['NOTE = Subject Series ' num2str(i) ' does not exist'])
            disp([A{i}])
            return
        else
        end
    end
    
    for i = 1:length(B)
        if rem(i,100)==0
            disp([' Validating file existence ' num2str(i)]);toc;
        end
        if exist(B{i},'file') == 0
            disp(['NOTE = Subject motion ' num2str(i) ' does not exist'])
            disp([B{i}])
            return
        else
        end
    end
end
%%Hardcode warning
%make an empty  matrix to write to:

all_overlaps_vec = zeros(91,(size(A,1)));
all_overlaps_vec_surf_only = zeros(91,(size(A,1)));

ncort_grey = 59412;

%reconcile two concs.
if subsetmotion ==1
    
    for i = 1: size(B,1)
        conc_sub_line = B{i};
        line_contents = strsplit(conc_sub_line,'/');
        sub_ses = char(line_contents(end));
        motion_sub{i,1} = char(sub_ses(1:23));
    end
    
    for i = 1: size(A,1)
        conc_sub_line = A{i};
        %overlap_sub{i,1} = char(conc_sub_line(88:110));
        line_contents = strsplit(conc_sub_line,'/');
        sub_ses = char(line_contents(end));
        overlap_sub{i,1} = char(sub_ses(1:23));
    end
    
    j=1;
    %find subjects from conc in csv:
    %try
        for i = 1: size(overlap_sub,1)
            motion_idx_cell{i,1} = find(contains(motion_sub,overlap_sub{i}));
            if size(motion_idx_cell{i,1},1) > 1 % find subjects that may have had repeated scans on the same day.
                rep_idx = motion_idx_cell{i,1};
                motion_idx(j,1) = rep_idx(1);j = j+1;
            else
                motion_idx(j,1) = motion_idx_cell{i,1};j = j+1;
            end
        end
    %catch
        %disp(i)
        %return
    %end
    motion_idx = unique(motion_idx);
end

for i =1:size(A)
    disp(num2str(i))
    cii = ciftiopen(A{i},wb_command);
    %cii = ciftiopen(overlap_dtseries_conc{i},wb_command);
    dtseries = cii.cdata;
    large_scalar_array_net = dtseries;
    large_scalar_array_net(:,4) = [];
    large_scalar_array_net(:,5) = [];
    
    nets = unique(nonzeros(large_scalar_array_net));
    
    overlap_grey_matrix = zeros(size(nets,1),size(nets,1));
    overlap_grey_matrix_surf_only = zeros(size(nets,1),size(nets,1));
    overlap_grey_matrix_surf_only_SA = zeros(size(nets,1),size(nets,1));
    [SA_greys,~,~] = surfaceareafromgreyordinates(C{i},D{i},1,A{i},'ADHD_overlapping_SA');
    %[all_areas_vec, network_surfarea, network_volume ] = surfaceareafromgreyordinates(Lmidthicknessfile,Rmidthicknessfile,output_only_greySA,dscalarwithassignments,outputname)
    %Get overlap of networks across by number of greyordinates.
    for n=nets(1):nets(end)
        nidx = find(dtseries(:,n)==n);
        for p=nets(1):nets(end)
            pidx = find(dtseries(:,p)==p);
            overlap_grey_matrix(n,p) = size(intersect(nidx,pidx),1);          
            %overlap_grey_matrix_SA(n,p) = sum(SA_greys(intersect(nidx,pidx)));
        end
    end
    
        %Get overlap of networks across surface area
    for n=nets(1):nets(end)
        nidx = find(dtseries(1:ncort_grey,n)==n);
        for p=nets(1):nets(end)
            pidx = find(dtseries(1:ncort_grey,p)==p);
            overlap_grey_matrix_surf_only(n,p) = size(intersect(nidx,pidx),1);
            overlap_grey_matrix_surf_only_SA(n,p) = sum(SA_greys(intersect(nidx,pidx)));          
        end
    end
    
    overlap_grey_matrix(:,4) =[];
    overlap_grey_matrix(4,:) =[];
    overlap_grey_matrix(:,5) =[];
    overlap_grey_matrix(5,:) =[];
    
    Ot=overlap_grey_matrix.';
    %m = (1:size(Ot,1)).' >= (1:size(Ot,2));
    m = triu(true(size(Ot)),1);
    Over_vec = Ot(m);
    
    overlap_grey_matrix_surf_only(:,4) =[];
    overlap_grey_matrix_surf_only(4,:) =[];
    overlap_grey_matrix_surf_only(:,5) =[];
    overlap_grey_matrix_surf_only(5,:) =[];
    
    overlap_grey_matrix_surf_only_SA(:,4) =[];
    overlap_grey_matrix_surf_only_SA(4,:) =[];
    overlap_grey_matrix_surf_only_SA(:,5) =[];
    overlap_grey_matrix_surf_only_SA(5,:) =[];   
    
    
    Ot=overlap_grey_matrix_surf_only.';
    %m = (1:size(Ot,1)).' >= (1:size(Ot,2));
    m = triu(true(size(Ot)),1);
    Over_vec_surf_only = Ot(m);
    
    Ot_SA=overlap_grey_matrix_surf_only_SA.';
    %m = (1:size(Ot,1)).' >= (1:size(Ot,2));
    m_SA = triu(true(size(Ot_SA)),1);
    Over_vec_surf_only_SA = Ot_SA(m);    
    
    %save matrix
    all_overlaps_square(:,:,i) =  overlap_grey_matrix;
    all_overlaps_square_surf(:,:,i) =  overlap_grey_matrix_surf_only;
    all_overlaps_square_surf_SA(:,:,i) = overlap_grey_matrix_surf_only_SA;
    
    %save vector
    all_overlaps_vec(:,i) = Over_vec;
    all_overlaps_vec_surf_only(:,i) = Over_vec_surf_only;
    all_overlaps_vec_surf_only_SA(:,i) = Over_vec_surf_only_SA; 
    
    % figure()
    % imagesc(overlap_grey_matrix); colormap jet
    % figure()
    % imagesc(triu(overlap_grey_matrix,1)); colormap jet
    %surface_overlaps(:,i) =Over_vec;
    disp('done')
end


figure()
subplot(plotcols,2,1)
imagesc(all_overlaps_vec); colormap jet; caxis(colormap_limits);
title('Subjects w/all frames of data.'); xlabel('subject');ylabel('net-netoverlap')



if subsetmotion ==1
    subjects_with_motion = B(motion_idx);
    TR_list = TRs(motion_idx);
    addpath('/home/exacloud/lustre1/fnl_lab/code/internal/analyses/compare_matrices/')
    [all_subjects_minutes,all_mean_FD] = check_twins_motion(subjects_with_motion,TR_list, 31,[],[]);
    
    subplot(plotcols,2,2)
    histogram(all_subjects_minutes, 20);
    title(['Subjects w/all available data. N=' num2str(size(all_subjects_minutes,1))]);xlabel('minutes');ylabel('count'); xlim([0 20]);
    
    for  j = 1:plotcols-1
        subplot(plotcols,2,(j*2)+1)
        all_cutoff_minutes_sub_idx = find(all_subjects_minutes>=minutes_cutoff(j));
        
        trimmed_subjects{j,1} = all_overlaps_vec(:,all_cutoff_minutes_sub_idx);
        imagesc(trimmed_subjects{j,1})
        
        title(['Subjects w/' num2str(minutes_cutoff(j)) '+ minutes of data.' ]); xlabel('subject');ylabel('net-netoverlap');caxis(colormap_limits);
        
        subplot(plotcols,2,(j*2)+2)
        histogram(all_subjects_minutes(all_cutoff_minutes_sub_idx), 20);
        title(['Subjects w/' num2str(minutes_cutoff(j)) '+ minutes of data. N=' num2str(size(all_cutoff_minutes_sub_idx,1))]);xlabel('minutes');ylabel('count'); xlim([0 20]);
        
        
        %save additional files just in case
        minutes_cutoff_all_overlaps_vec{j,1} = all_overlaps_vec(:,all_cutoff_minutes_sub_idx);
        minutes_cutoff_all_overlaps_vec_surf_only{j,1} = all_overlaps_vec_surf_only(:,all_cutoff_minutes_sub_idx);
        minutes_cutoff_subject_list{j,1} = subjects_with_motion(all_cutoff_minutes_sub_idx);
        minutes_cutoff_all_overlaps_square{j,1} = all_overlaps_square(:,:,all_cutoff_minutes_sub_idx);
        minutes_cutoff_all_overlaps_square_surf{j,1} = all_overlaps_square_surf(:,:,all_cutoff_minutes_sub_idx);
        cutoff_indices{j,1} = all_cutoff_minutes_sub_idx;
        minutes_cutoff_all_overlaps_vec_surf_only_SA{j,1} = all_overlaps_vec_surf_only_SA(:,all_cutoff_minutes_sub_idx);
        minutes_cutoff_all_overlaps_square_surf_SA{j,1} = all_overlaps_square_surf_SA(:,:,all_cutoff_minutes_sub_idx);
    end
end
save('all_trio_and_prisma_TM_cleaned_numgreyordinates_in_overlaps_w_surfarea.mat','all_overlaps_vec_surf_only_SA','all_overlaps_vec','all_overlaps_vec_surf_only','all_overlaps_square','all_overlaps_square_surf','minutes_cutoff_subject_list','minutes_cutoff','minutes_cutoff_all_overlaps_vec','minutes_cutoff_all_overlaps_vec_surf_only','minutes_cutoff_all_overlaps_square','minutes_cutoff_all_overlaps_square_surf','minutes_cutoff_all_overlaps_vec_surf_only_SA','minutes_cutoff_all_overlaps_square_surf_SA','cutoff_indices')


end