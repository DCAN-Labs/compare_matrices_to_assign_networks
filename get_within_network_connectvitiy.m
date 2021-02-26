function [ output_args ] = get_within_network_connectvitiy
%This function calculates the within-netwwork connectivity for
%both the template matching and gordon pconns.

%   Detailed explanation goes here

gordon =0;
both_methods =0;

% this_code = which('simple_pconn_average');
% [code_dir,~] = fileparts(this_code);
% support_folder=[code_dir '/support_files']; %find support files in the code directory.
% addpath(genpath(support_folder));
% settings=settings_comparematrices;%
% np=size(settings.path,2);
%
% disp('Attempting to add neccesaary paths and functions.')
% warning('off') %supress addpath warnings to nonfolders.
% for i=2:np
%     addpath(genpath(settings.path{i}));
% end
% rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
% rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
% addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
% warning('on')
% wb_command=settings.path_wb_c; %path to wb_command
wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
if gordon ==1
    %cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/projects/ABCD/avg_pconn_maker/group2_10min_mean.pconn.nii',wb_command);
    %parcel_file =load('/home/exacloud/lustre1/fnl_lab/projects/Polyvertexscore/HumanGordon_parcel.mat');
    cii1 = ciftiopen('C:\Users\hermosir\Desktop\TM_pconns\group1_10min_mean.pconn.nii',wb_command);
    cii2 = ciftiopen('C:\Users\hermosir\Desktop\TM_pconns\group2_10min_mean.pconn.nii',wb_command);
    parcel_file =load('C:\Users\hermosir\Desktop\TM_pconns\HumanGordon_parcel.mat');
            parcel = parcel_file.parcel; % for TM

    parcel1 = parcel_file.parcel; % for gordon
    pconn1 = cii1.cdata;
    pconn2 = cii2.cdata;
    pconn1 = tanh(pconn1);
    pconn2 = tanh(pconn2);
    
    netsprequeeze=struct2cell(parcel1); % convert struct to cell to get parcel names
    
    net_names = squeeze(netsprequeeze);
    
    %cii = ciftiopen('/home/exacloud/lustre1/fnl_lab/projects/ABCD_net_template_matching/pconns/ABCD_threshold_template_maps_grp2_average_pearson.pconn.nii',wb_command);
    %parcel_file =load('/home/exacloud/lustre1/fnl_lab/projects/Polyvertexscore/parcel_probability_map.mat');
else
    cii3 = ciftiopen('C:\Users\hermosir\Desktop\TM_pconns\ABCD_threshold_template_maps_grp1_average.pconn.nii',wb_command);
    cii4 = ciftiopen('C:\Users\hermosir\Desktop\TM_pconns\ABCD_threshold_template_maps_grp2_average.pconn.nii',wb_command);
    parcel_file =load('C:\Users\hermosir\Desktop\TM_pconns\parcel_probability_map.mat');
    
    
    if both_methods ==1
        parcel2 = parcel_file.parcel; % for TM
        netsprequeeze=struct2cell(parcel2); % convert struct to cell to get parcel names
        net_names = squeeze(netsprequeeze);
    else
        
        parcel = parcel_file.parcel; % for TM
        netsprequeeze=struct2cell(parcel); % convert struct to cell to get parcel names
        net_names = squeeze(netsprequeeze);
    end
end

if both_methods ==1
    pconn3 = cii3.cdata;
    pconn4 = cii4.cdata;
else
    

    
    if gordon ==0
               	pconn1 = cii3.cdata;
        pconn2 = cii4.cdata;
        pconn1 = tanh(pconn1);
        pconn2 = tanh(pconn2);
    else
       	pconn1 = cii1.cdata;
        pconn2 = cii2.cdata;
    end
end




if both_methods ==1
    pconns{:,:,1} = pconn1;
    pconns{:,:,2} = pconn2;
    pconns{:,:,3} = pconn3;
    pconns{:,:,4} = pconn4;
else
    pconns{:,:,1} = pconn1;
    pconns{:,:,2} = pconn2;
end


%triupconn = triu(pconn,1);
concat_connection_origins = [];
concat_all_connections = [];
connection_origins = [];
all_connections = [];


%net_names_both = [net_names_gordon;net_names_TM];


for p = 1:size(pconns,3)
    pconn = pconns{:,:,p};
    if both_methods ==1
        if  p>2
            parcel = parcel2;
        else
            parcel = parcel1;
        end
    else
    end
    for i = 1:size(parcel,2)
        this_parcels_struct = parcel(i);
        mini_conn = pconn(this_parcels_struct.ix,this_parcels_struct.ix); % get only the connections for this network.
        triu_mini_conn_vec = mini_conn(triu(true(size(mini_conn,1)),1)); % unwrap upper triangle into a vector.
        connections_by_network{i,1} = triu_mini_conn_vec;
        %netname = parcel.shortname(i);
        %netname = char(net_names(2,i)); % line 2 is where the short names live.
        if both_methods ==1
            if p <3
                netname = char(net_names_gordon(2,i)); % line 2 is where the short names live.
                repeatednamesforthisnetwork=repmat(netname,size(triu_mini_conn_vec,1),1);
            else
                netname = char(net_names_TM(2,i)); % line 2 is where the short names live.
                repeatednamesforthisnetwork=repmat(netname,size(triu_mini_conn_vec,1),1);
            end
        else
            netname = char(net_names(2,i)); % line 2 is where the short names live.
            repeatednamesforthisnetwork=repmat(netname,size(triu_mini_conn_vec,1),1);
            
        end
        if i==1
            connection_origins = repeatednamesforthisnetwork;
            all_connections = triu_mini_conn_vec;
        else
            connection_origins = strvcat(connection_origins, repeatednamesforthisnetwork);
            all_connections = [all_connections; triu_mini_conn_vec];
        end
        
        average_within_net_connextivity(i,1) = mean(triu_mini_conn_vec);
        std_within_net_connextivity(i,1) = std(triu_mini_conn_vec);
        
        
    end
    if p ==1
        concat_connection_origins = connection_origins;
        concat_all_connections = all_connections;
    else
        concat_connection_origins = [concat_connection_origins; connection_origins];
        concat_all_connections = [concat_all_connections; all_connections];
    end
end

connect_cells=cellstr(concat_connection_origins);
%each was run seperately
%group1_Gordon_mean = average_within_net_connextivity;
%group2_Gordon_mean = average_within_net_connextivity;
%group1_TM_mean = average_within_net_connextivity;
%group2_TM_mean = average_within_net_connextivity;
%group1_Gordon_std = std_within_net_connextivity;
%group2_Gordon_std = std_within_net_connextivity;
%group1_TM_std = std_within_net_connextivity;
%group2_TM_std = std_within_net_connextivity;

%if grouping_var ==1
grouping_var = [repmat(1,size(all_connections,1),1);repmat(2,size(all_connections,1),1)];
%end

disp('Wait here.')
F=figure();
net_colors = cell2mat(net_names(5,1:end)');
net_colors = cell2mat(net_names(5,1:end)');
if gordon ==0
    net_colors = net_colors (:,1:3);
end

if both_methods ==1
for c= 1: size(net_colors,1)*2
    double_net_colors(c,:) = net_colors(round(c/2),:);
end
end
%boxplot(all_connections, connection_origins)
G = boxplot(concat_all_connections, {concat_connection_origins,grouping_var},'Symbol','o','Colors',net_colors,'jitter',0.5);
%G = boxplot(concat_all_connections, concat_connection_origins,'Symbol','o','Colors',net_colors,'jitter',0.5);

MarkerFaceColor(net_colors)
save('TMandGordon_within_network_stats.mat','group1_TM','group2_TM','group1_Gordon','group2_Gordon')
net255=net_colors*255;



%example
data = rand(20,24);
month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];

end

