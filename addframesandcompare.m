function addframesandcompare(input_series,mask,num_distributed_command)

%This function generates connectivity matrices and 
%Inputs


%% Load settings
support_folder=[pwd '/support_files'];
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

file_split = strsplit(file_dir,'/');
file_name = char(file_split(end));

for i=1:np
    addpath(genpath(settings.path{i}));
end

wb_command=settings.path_wb_c; %path to wb_command


conc = strsplit(input_series, '.');
conc = char(conc(end));

if strcmp('conc',conc) == 1
    A = importdata(input_series);
    for i = 1:length(A) %check for subjects
        if exist(A{i}) == 0
            disp(['Subject series ' num2str(i) ' does not exist'])
            %return
        else
        end
    end
else
    A = {input_series};
    for i = 1:length(A) %check for subjects
        if exist(A{i}) == 0
            disp(['Subject series does not exist'])
            %return
        else
        end
    end 
end

if strcmp('conc',conc) == 1
    B = importdata(mask);
    for i = 1:length(B) %check for subjects
        if exist(B{i}) == 0
            disp(['Subject mask ' num2str(i) ' does not exist'])
            %return
        else
        end
    end
else
    B = {mask};
    for i = 1:length(B) %check for subjects
        if exist(B{i}) == 0
            disp(['Subject mask does not exist'])
            %return
        else
        end
        mask_vectorinv = load(B{1});
        mask_vector = ~logical(mask_vectorinv');
    end 
end

    
    disp(A{i});
    cii=ciftiopen(A{i},wb_command);
    newcii=cii;
    pt_series=single(newcii.cdata);
    N_rois=size(pt_series,1);% count the total rois in this parcellation
    all_framesmatrix = single(zeros(N_rois,N_rois));
    
    tic
    %build all frames matrix for comparison
    for r=1:size(pt_series,1)
        current_row = pt_series(r,:);
        for c=1:size(pt_series,1)
           current_column = pt_series(c,:);
            all_framesmatrix(r,c) = corr(current_row(mask_vector)',(current_column(mask_vector))');            
        end
    end
    toc
    
    iterative_matrix= single(zeros(N_rois,N_rois,(size(pt_series,2)-1))); % create 3D matrix that is ROI X ROI X frames included
    
    
    if num_distributed_command > 0%disribute srun commands across nodes
        command_indices = 1:num_distributed_command;
        dir1 = [pwd '/' filename];
        if exist(dir1)==0
            cmd = ['mkdir ' dir1];
            system(cmd);
        end
        
        %The code will split the subject matrix into a series of columns of width.  
        %The code will the attempt to correlate each column in the split subject matrix with the each of the columns assinged in the all frames matrix.
        
        for i = 1:num_distributed_command
        cmd = ['srun --time=36:00:00 --mem=64Gb -e ' filename i '.err -o ' filename i '.out run_distributed_correlation.sh ' input_vectors '_i' all_frames_matrix_path ' output_pieces_i.m' ]; 
        system(cmd);
        end
        
    else %handle all this inside matlab.
        
        %build matrices from smaller masks
        for m=2:size(pt_series,2)
            current_mask_vector = mask_vector(1:m);
            for r=1:size(pt_series,1)
                current_row = pt_series(r,:);
                for c=1:size(pt_series,1)
                    current_column = pt_series(c,:);
                    iterative_matrix(r,c,m) = corr(current_row(current_mask_vector)',(current_column(current_mask_vector))');
                    %R=corr(pt_series_final(:,mask)',dt_series(:,mask)');
                end
            end
        end
    end
    
    
    for i=2:size(pt_series,2)
        
    end
    
    A=size(inputseries,1);
    mat_size = 352; %hardcode - number of parcellations
    lowern=mat_size-1; %used for lower calculation of triangle
    Gausse_num=((lowern*(lowern+1))/2);
    subjectmat=zeros(Gausse_num,A);
    
end
