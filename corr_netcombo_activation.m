function corr_netcombo_activation(probabilistic_network_conc,activation_dscalar,wb_command,maximum_combination_of_nets,outputfolder,use_negative_nets, abs_activation, corr_type, create_max_combo_cifti, output_name ,run_locally)

% R. Hermosillo
% v1.0 - 07/28/2021

% This function tries to find a combination of probabilistic dscalars that
% has maximal correlation with task activation map.  The general theory being
% that activation patterns are combination of network activations. It does
% this by calculating the addative sum of all possible combinations of
% networks, then calculates a ranked correlation (Kendall's Tau) with the
% activation map.  Kendall's Tau is used (rather than Pearson) because all probabailistic maps
% are positive and the activation map values are postive and negative.

%Inputs are:
% probabilistic_network_conc = Char. a confile (list of networks) IMPORTANT NOTE:
% this conc MUST be in the same order as the Network number ("net_list" variable below). (i.e. DMN
% first, Vis second, etc.).  If not, the combinations of "network 1" and
% "network 7" will actually correspond with the addative sum of different
% networks.
 
%activation_dscalar = Char. a dscalar file of activation (This should be the
% same (i.e. 91282 greyordinates in length as the probabilistic dscalars).

% wb_command = Char. path to your install version of wb_command.
% maximum_combination_of_nets = Interger. Limit the maximum number of combinations.

% use_negative_nets = While this function tries to find the additive
% combination of networks that correspond with activation, it's also
% likely that some networks are actively supressed.  If you set this to 1,
% the function will also try every combination of networks, with and some
% subtracted. 0 assumes networks are additive.

% abs_activation = take the absolute value of the activation map.
% corr_type = Provide 'Kendell' or 'Pearson'

% creat_max_combo_cifti = build a cifti that is based off of the maximum
% correlation combination. This usese the activation cifti as the template.

% output_name = The name of you output file.  THis will be appended with
% the selected parameters.

% run_locally = Set to 0 for running on MSI.  Set to 1 for Robert's computer to debug.  
% 

%Additional packages are required to handle ciftis:
%MSCcodebase-master: https://github.com/MidnightScanClub/MSCcodebase
%HCP packages: https://github.com/Washington-University/HCPpipelines/

%Example call:
%corr_netcombo_activation('C:\Users\...probability_maps_in_order.conc','C:\Users\...\activation_map.dscalar.nii','C:\Users\...\wb_command', 5, 'C:\Users\...someoutputdirectory', 1, 0, 'Pearson', 1, 'myfilename',0);

% Step 0: check is potentially numeric variables are unintentionally imported from the command line as strings:
if isnumeric(maximum_combination_of_nets) ==1
else
    maximum_combination_of_nets = str2num(maximum_combination_of_nets);
end
if isnumeric(use_negative_nets) ==1
else
    use_negative_nets = str2num(use_negative_nets);
end
if isnumeric(abs_activation) ==1
else
    abs_activation = str2num(abs_activation);
end

if isnumeric(create_max_combo_cifti) ==1
else
    create_max_combo_cifti=str2num(create_max_combo_cifti);
end
if isnumeric(run_locally)==1
else
    run_locally=str2num(run_locally);
end

%% Step 1: Add dependency paths

%add cifti paths
if run_locally ==1
    %Some hardcodes:
    %wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
else    
    this_code = which('template_matching_RH');
    [code_dir,~] = fileparts(this_code);
    support_folder=[code_dir '/support_files']; %find support files in the code directory.
    addpath(genpath(support_folder));
    settings=settings_comparematrices;%
    np=size(settings.path,2);
    disp('Attempting to add neccesaary paths and functions.')
    warning('off') %supress addpath warnings to nonfolders.
    %for i=2:np
    for i=1:np
        
        addpath(genpath(settings.path{i}));
    end
    warning('on')
    wb_command=settings.path_wb_c; %path to wb_command
    load([code_dir filesep 'support_files' filesep 'PowerColorMap_wzero.mat'],'mycmap'); %load the power color map
end

%% Step 2: Load activation data
% 
net_list = [1 2 3 5 7 8 9 10 11 12 13 14 15 16]; % hardcoded network assingments.
all_labels = {'DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};


activation_dscalar_file = ciftiopen(activation_dscalar,wb_command);
activation_dscalar_cifti = activation_dscalar_file.cdata;

if abs_activation ==1
    activation_dscalar_cifti = abs(activation_dscalar_cifti);
    abs_act = 'ON'; % for saved file name
else
    abs_act = 'OFF';
end

%reference_cifti = reference_cifti_all(1:59412,1); %just use surface labels
net_file_list = importdata(probabilistic_network_conc);
num_nets = size(net_list,1);

%% Step 3: Build list of combinations of networks
D=[]; k=1;

%for i = net_list(:)'
for i = 1:maximum_combination_of_nets
    C = nchoosek(net_list',k); % return all possible combinations of select k nets from all networks.
    % put the vector of C into a matrix of zeros so that it can be concatenated onto the full list.
    Czeromat = zeros(size(C,1),size(net_list,2),1);
    Czeromat(:,1:size(C,2)) = C;
    D = [D; Czeromat];
    k=k+1;
end

disp(['Based on the maximum combinations you have selected: ' num2str(maximum_combination_of_nets) ', there are '  num2str(size(D,1)) ' combinations.' ])

%build negative combinations
if use_negative_nets ==1
    for N = 1:maximum_combination_of_nets
    binary_combinations_all{1,N} = dec2bin(0:2^N-1)-'0'; % get all binary combinations of max values.
    zeroindxs = binary_combinations_all{1,N} ==0;
    binary_combinations_all{1,N}(zeroindxs) = -1; % set 0 values to neagtive 1.  So networks will be multiplied by 1 or negative 1.
    end
    negnets= 'ON';
else
    negnets= 'OFF';
end

%% Step 4: Load each probability cifti then build network probabiliity matrix
n = 1;m=1;
%for i = 1:size(net_file_list,1)
for i = 1:max(net_list)
    if i ==4 || i==6 % HARDCODE Warning: In template matching, these numbered networks don't exist.
        net_prob(:,n) = nan(size(net_prob,1),1); %probably 91282 
    else
        ciifile = ciftiopen(net_file_list{m},wb_command);
        net_prob(:,n) = ciifile.cdata;
        m=m+1;
    end
    n=n+1;
end

%% Step 5: Start correlation of activation with each combination of network.

tic
if use_negative_nets ==1
    summation_tau_vec = zeros(size(D,1),size(binary_combinations_all{1,end},1)); % preallocate for speed.
else
    summation_tau_vec = zeros(size(D,1),1); % preallocate for speed.
end

for i = 1:size(D,1) % use 3472 for maximum 5 nets.
    
    if rem(i,100) ==0 %give updates every 100 calculations.
        disp(i);toc;
    end
    
    nets_combo =net_prob(:,nonzeros(D(i,:)));
    
    if use_negative_nets ==1
        current_number_of_nets = size(nonzeros(D(i,:)),1);
        this_number_of_binary_combinations =binary_combinations_all{1,current_number_of_nets};
        for N = 1:size(this_number_of_binary_combinations,1)
            nets_combo_flipping = nets_combo.*(this_number_of_binary_combinations(N,:));
            netssum = sum(nets_combo_flipping,2);
            switch corr_type
                case 'Kendall'
            [summation_tau_vec(i,N),~] = corr(activation_dscalar_cifti,netssum,'Type','Kendall');
                case 'Pearson'
                    [summation_tau_vec(i,N),~] = corr(activation_dscalar_cifti,netssum);
                otherwise
                    error('Correlation type not supported. Select Kendall or Pearson.')
            end
        end
    else
        netssum = sum(nets_combo,2);
        %kendall_tau_vec(i,1) = kendall_tau(activation_dscalar_cifti,netssum);
        %%don't use external kendall_tau function. it's  much slower.
        switch corr_type
            case 'Kendall'
                [summation_tau_vec(i,1),~] = corr(activation_dscalar_cifti,netssum,'Type','Kendall');
            case 'Pearson'
                [summation_tau_vec(i,1),~] = corr(activation_dscalar_cifti,netssum);
            otherwise
                error('Correlation type not supported. Select Kendall or Pearson.')
        end
    end
end

%% Step 6 Find max and Report combination
if use_negative_nets ==1
    maxtau = max(max(summation_tau_vec));
    [max_tauindx,binary_combinationidx]=find(summation_tau_vec==maxtau);
    
else
    [maxtau,max_tauindx]= max(summation_tau_vec);
    
end
max_combo = nonzeros(D(max_tauindx,:))';
combo_names = all_labels(max_combo);

disp(['Combination of networks that yeilds the highest rank corelation is number: ' num2str(max_combo)  ' = ' num2str(maxtau)]);
disp(['...using ' num2str(size(nonzeros(D(max_tauindx,:)),1)) ' network(s)' ])
if use_negative_nets ==1
    disp(['The positive-negative ordering of that combination is : ' num2str(binary_combinationidx)]);
end

disp('Combination of networks that yeilds the highest rank corelation is : ' );
disp(combo_names)
disp('Saving correlation output and parameters to file...')

for j=1:size(combo_names,2)
    if j==1
        combo_char = char(combo_names{j});
    else
        combo_char_next = char(combo_names{j});
        combo_char=[combo_char combo_char_next];
    end
end

%% Step 7 Save data (.mat file with paramters)
if use_negative_nets ==0
    save([outputfolder filesep output_name '_corrtype' corr_type '_maxnets' num2str(maximum_combination_of_nets) '_absactivation' abs_act '_negnets' negnets '_max_probnet_combination_' combo_char '.mat'], 'D', 'summation_tau_vec', 'probabilistic_network_conc','activation_dscalar','wb_command','maximum_combination_of_nets','combo_names','outputfolder', 'use_negative_nets', 'abs_activation', 'corr_type', 'create_max_combo_cifti', 'run_locally');
else
    save([outputfolder filesep output_name '_corrtype' corr_type '_maxnets' num2str(maximum_combination_of_nets) '_absactivation' abs_act '_negnets' negnets '_max_probnet_combination_' combo_char '.mat'], 'D', 'summation_tau_vec', 'probabilistic_network_conc','activation_dscalar','wb_command','maximum_combination_of_nets','combo_names','outputfolder', 'use_negative_nets', 'abs_activation', 'corr_type', 'create_max_combo_cifti', 'run_locally','binary_combinationidx');
end

%% Step 8 (Optional): Save combination dscalar.
if create_max_combo_cifti ==1
    disp('Saving combination file...')
    if use_negative_nets ==0
        max_nets_combo =net_prob(:,max_combo);
        netssum = sum(max_nets_combo,2);
    else
        % first get the combination that was found to be the max
        max_nets_combo = net_prob(:,nonzeros(D(max_tauindx,:)));
        %next use the size of the combination to get reference the binary
        %sets
        max_binary_flip_set = binary_combinations_all{1,size(combo_names,2)};
        % Use the max_binary index to find the -1/1 combination in the set.
        max_binary_flip = max_binary_flip_set(binary_combinationidx,:);
        
        %invert combinations
        nets_combo_flipping = max_nets_combo.*(max_binary_flip);
        %sum
        netssum = sum(nets_combo_flipping,2);
    end
    
    ciifile.cdata = netssum;
    ciftisave(ciifile,[outputfolder filesep output_name '_corrtype' corr_type '_maxnets' num2str(maximum_combination_of_nets) '_absactivation' abs_act 'negnets' negnets '_max_probnet_combination_' combo_char '.dscalar.nii'],wb_command)
end

figure()
if use_negative_nets ==0
    set(gcf,'color','w');
    %imagesc(D(1:3472,1:5)');
    imagesc(D(:,1:maximum_combination_of_nets)');
    colormap(mycmap);
    xlabel('Combination')
    ylabel('Number of networks in combination')
    hold on
    yyaxis right
    p = plot(summation_tau_vec,'LineWidth',2,'Color','k','LineStyle','-');
    ylabel('Correlation')
    hold on
    plot(max_tauindx,maxtau,'wo','MarkerSize',7,'MarkerFaceColor','w')
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    title(output_name)
    
    % Make sorted polot
    figure;
    set(gcf,'color','w');
    [corr_sorted, corr_sorted_idx] = sort(summation_tau_vec,'descend');
    D_sorted= D(corr_sorted_idx,:);
    imagesc(D_sorted(:,1:maximum_combination_of_nets)');
    colormap(mycmap);
    xlabel('Combination')
    ylabel('Number of networks in combination')
    hold on
    yyaxis right
    p = plot(corr_sorted,'LineWidth',2,'Color','k','LineStyle','-');
    ylabel('Correlation')
    hold on
    plot(max_tauindx,maxtau,'wo','MarkerSize',7,'MarkerFaceColor','w')
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    title(output_name)
    
else
    disp('Plotting data for negative nets not yet implemented.')
end
%disp(['Saving image: ' opts.saveFolder outputfilename])
%print([opts.saveFolder outputfilename], '-dpng', '-r600')

disp('Done running correlation script.')