function corr_netcombo_activation(probabilistic_network_conc,activation_dscalar,wb_command,maximum_combination_of_nets,outputfolder)

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
% activation_dscalar = Char. a dscalar file of activation (This should be the
% same (i.e. 91282 greyordinates in length as the probabilistic dscalars).
% wb_command = Char. path to your install version of wb_command.
% maximum_combination_of_nets = Interger. Limit the maximum number of combinations.

%Additional packages are required to handle ciftis:
%MSCcodebase-master: https://github.com/MidnightScanClub/MSCcodebase
%HCP packages: https://github.com/Washington-University/HCPpipelines/

%Example call:
%corr_netcombo_activation('C:\Users\...probability_maps_in_order.conc','C:\Users\...\activation_map.dscalar.nii','C:\Users\...\wb_command',5);


%% Step 1: Add dependency paths

%add cifti paths
%wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')

%load power colors
load('C:\Users\hermosir\Documents\repos\support_folder\Jet_wzerowhite_colormap.mat','mymap')

%% Step 2: Load activation data

net_list = [1 2 3 5 7 8 9 10 11 12 13 14 15 16]; % hardcoded network assingments.
all_labels = {'DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};


activation_dscalar_file = ciftiopen(activation_dscalar,wb_command);
activation_dscalar_cifti = activation_dscalar_file.cdata;
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

%% Step 5: Start ranked correlation of activation with each combination of network.
tic
for i = 1:size(D,1) % use 3472 for maximum 5 nets.
    if rem(i,100) ==0 %give updates every 100 calculations.
        disp(i);toc;
    end
    nets_combo =net_prob(:,nonzeros(D(i,:)));
    netssum = sum(nets_combo,2);
    %kendall_tau_vec(i,1) = kendall_tau(activation_dscalar_cifti,netssum);
    %%don't use kendall_tau function. it's  much sower
    [tau_vec(i,1),~] = corr(activation_dscalar_cifti,netssum,'Type','Kendall');
end


[maxtau,max_tauindx]= max(tau_vec);

max_combo = nonzeros(D(max_tauindx,:))';

combo_names = all_labels(max_combo);
disp(["Combination of networks that yeilds the highest rank corelation is : " num2str(max_combo)  ' = ' num2str(maxtau)]);
disp("Combination of networks that yeilds the highest rank corelation is : ");
disp(combo_names)
%% Step 6: Save combination dscalar.
if creat_max_combo_cifti ==1
    disp('Saving combination file.')
    max_nets_combo =net_prob(:,nonzeros(D(max_tauindx,:)));
    netssum = sum(nets_combo,2);
    ciifile.cdata = netssum;
    ciftisave(ciifile,[outputfolder filesep 'max_net_citfti_combination' combo_names '.dscalar.nii'],wb_command)
end

% figure()
% imagesc(jaccard_matrix_by_subject_by_network')
% set(gcf,'color','w');
% xlabel('subject')
% ylabel('TM networks')