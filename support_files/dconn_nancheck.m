function cifti_nancheck(cifti_conc_file)

%R.hermosillo This function checks for nans in your cifti
%4/10/2019 R. Hermosillo added an option to pass in a conc file of ciftis.

conc = strsplit(cifti_conc_file, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    ciftis = importdata(cifti_conc_file);
else
    ciftis = {cifti_conc_file};
end

for i = 1:length(ciftis)
    if exist(ciftis{i}) == 0
        disp(['NOTE = Subject Series ' num2str(i) ' does not exist'])
        return
    else
    end
end
disp('All series files exist continuing ...')
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/Matlab_CIFTI'))
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/gifti-1.6'))

%path_wb_c='/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-9253ac2/bin_rh_linux64/wb_command';
path_wb_c='/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-1.2.3-HCP/bin_rh_linux64/wb_command';
cifti_type = strsplit(ciftis{1}, '.');
cifti_exten = char(cifti_type(end-1));

if strcmp('dtseries',cifti_exten) == 1
    
elseif strcmp('dconn',cifti_exten) == 1
    
elseif strcmp('dscalar',cifti_exten) == 1
    
else
    disp('filetype not supported by not checker: ')
end

k = 1;
for i = 1:length(ciftis)
    newcii = ciftiopen(ciftis{i}, path_wb_c);
    
    cifti=single(newcii.cdata);
	%[m, n] = size(cifti);
    clear newcii
    disp(i)
    %disp(row);
    %disp(col);
    switch cifti_exten
        
        case 'dtseries' % for use with dtseries. % goes through every time point  
            [row, col] = find(isnan(cifti));
            if isempty(row) ==1
            else
                disp(['This subject ' num2str(i) 'has Nans'])
                disp(['File with nans is : ' ciftis{i} ])
            end

	%for j = 1:n
        %if sum(isnan(cifti(:,j))) > 0
        %    disp(['This subject ' num2str(i) ' has Nans'])
        %    disp(['File with nans is : ' ciftis{i} ])
        %    disp(['number of greyordinates with nans = ' num2str(sum(isnan(cifti(:,j))))])
        %    badsubidx(k,1) = i; k = k+1;
        %else
        %end
	%end

        otherwise % for use with dconn or dscalar.
            if sum(isnan(cifti(:,1))) > 0
                disp(['This subject ' num2str(i) 'has Nans'])
                disp(['File with nans is : ' ciftis{i} ])
                disp(['number of greyordinates with nans =' sum(isnan(cifti(:,1)))])  
            else
            end
    end
end
