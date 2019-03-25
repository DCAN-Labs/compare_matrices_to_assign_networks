
%load conc file of test dscalars.
%A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztemplate_match_intervals/MSC_sofar.conc');
A = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Ztempl_match_and_infomap_intervals/infomap_all_subjects1run.conc');
%load confile of reference dscalars.
B = importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/template_matching_none_minutes_half2.conc');
%B= importdata('/home/exacloud/lustre1/fnl_lab/projects/midnight_scan_club/template_matching/Zscored_templatematching/infomap_all_subjects1run.conc');
num_subjects = 10;
num_time_intervals = 8;
none_minuteslimit = 0;
num_reps = 1;

%% Adding paths for this function
this_code = which('plot_template_matching_reps');
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
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')

wb_command=settings.path_wb_c; %path to wb_command

%% Do checks
for i = 1:length(A)
    if exist(A{i}, 'file') == 0
        disp(['NOTE: Subject Series ' num2str(i) ' does not exist']);
        return
    else
    end
end
disp('All test scalars files exist continuing ...')

for i = 1:length(B)
    if exist(B{i}, 'file') == 0
        disp(['NOTE: Subject Series ' num2str(i) ' does not exist']);
        return
    else
    end
end
disp('All reference scalars files exist continuing ...')

if none_minuteslimit == 1
    expected_num_scalars = (num_time_intervals * num_reps * num_subjects) + num_subjects;
else
    expected_num_scalars = (num_time_intervals * num_reps * num_subjects);
end

if size(A,1) == expected_num_scalars
    disp('Number of subjects and reps matches conc file')
else
    disp('Number of subjects/reps is calculation does not match conc file.')
    return
end

if size(B,1) ~= num_subjects
    disp('Number of subjects in the conc file does not match the number of subjects')
    return
else
    disp('Number of subjects in the conc file does matches the number of subjects....continuing')
end
 m = 1; % counter for line in conc file.
%% Step 4: load template networks from 2nd half of data:
for i = 1:num_subjects
    reference_cii = ciftiopen(B{i},wb_command);
    other_half_networks = reference_cii.cdata;
    disp(['Reference scalar = ' B{i}])
    for j = 1:num_time_intervals+1
        if j <= num_time_intervals %check if the interval is the last one (i.e. "none minutes limit" for which there will be no reps.)
            for k = 1:num_reps
                test_cii = ciftiopen(A{m},wb_command);
                new_subject_labels = test_cii.cdata;
                disp(['Test scalar = ' A{m}]);
                disp('Caluating mutual information between scalars')
                %% Step 5: Calculate mutual information
                muI(i,j,k) = MutualInformation(new_subject_labels, other_half_networks); %Mutual information
                [VIn(i,j,k), MIn(i,j,k)] = partition_distance(new_subject_labels, other_half_networks); %Normalized variation of information ([p, q] matrix), Normalized mutual information ([p, q] matrix
            m =m+1;
            end
        else % if it is the minutes limit, dont go through reps.
            m =m+1;
        end
    end
    disp(num2str(i))
end
rep1 = squeeze(muI(:,:,1)');
plot(minutes, rep1);
C = permute(muI,[1,2,3]);
D = reshape(C,[],size(muI,2),1);
%groups = [[1:10]'; [1:10]'; [1:10]';];
groups = num2str([[1:10]'; [1:10]'; [1:10]';]);
subjects = cellstr(groups);
Dgroups = repmat([1:10]',1,8,3);
rep1b = squeeze(muI(:,:,1)); rep2b = squeeze(muI(:,:,2)); rep3b = squeeze(muI(:,:,3));E = [rep1b;rep2b;rep3b;];
disp('Done running plot_templatematching_reps');