function plotdconn(dconn_cifti_path,net_assigns,downsample_dconn,DS_factor,apply_Zscore_dconn,image_name,plot2dconns,dconn_cifti_path2,use_nets1,net_assigns2_file,caxis_scale,Pos_neg_colormap,exclude_zero_networks,save_non_zeromatrix,only_save_abs_diff_dconn)

%R.Hermosillo 08/20/2019

%This code make a picture of the dconn matrix, sorted by networks, the
%dconn is down sampled (DS_factor)to reduce the number of pixels in the image.

%inputs are:
% 1)path to a dconn
% 2)path to a dscalar with asssignments (ordinal).
% 3)downsample - set to 1 to down sample. set to 0 to plot the full dconn.
% 4)the downsampling factor. set to 1 to use the full data set.
% 5)Zscore dconn to to tre is you want to apply a within-region Zscore
%transformaiton of the data.
%6) image name
%7) if you want ot plot 2 dconns, set to 1.
%8) path to dconn2
%9) path to dscalar 2.
%11) use_nets1  = use the same networks to sort both dconns.
%12) net_assigns2_file = if the "use nets1" is 0, the supply the network assigns file.
%13) caxis_scale = the color axis for the dconns.
%14) Pos_neg_colormap  = use if you want to use a 0 equals white colormap,
%set to 1. 0 uses color map jet.
%15) exclude_zero_networks = to 1 if your dscalar has zeros and you want to
%exclude those grarordinates.

%net_order = [12 9 5 1 3 14 15 16 8 10 11 13 7 2];

%net_order = [10 7 4 1 3 12 13 14 6 8 9 11 5 2];
%Zscore_dconn=1; % set to true if you want to Zscore your dconn.
%downsample_dconn  =1;
%DS_factor = 50;

%% add dependencies
%parameters:
%wb_command='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command';
wb_command='/home/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command';
addpath(genpath('/panfs/jay/groups/6/faird/shared/code/external/utilities/gifti-1.6'))
addpath(genpath('/panfs/jay/groups/6/faird/shared/code/internal/utilities/Matlab_CIFTI'))
%load('/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/PowerColorMap_wzero.mat');
%load('parcel_probability_map.mat','parcel'); [~,index] = sortrows([parcel.power_val].'); parcel = parcel(index); clear index

if downsample_dconn == 1
    if isnumeric(DS_factor) ==1
        if DS_factor <=0
            error('DS_factor must be a positive real number.')
        end
    else
        error('DS_factor must be numeric')
    end
end

%load assignments

if strcmp(net_assigns(end-3:end),'.csv')
    assigns = table2array(readtable(net_assigns));
elseif strcmp(net_assigns(end-3:end),'.txt')
    assigns = table2array(readtable(net_assigns)); 
elseif strcmp(net_assigns(end-10:end),'pscalar.nii')
    net_assigns = ciftiopen(net_assigns,wb_command);
    assigns = net_assigns.cdata;
elseif strcmp(net_assigns(end-10:end),'dscalar.nii')
    net_assigns = ciftiopen(net_assigns,wb_command);
    assigns = net_assigns.cdata;
elseif isnumeric(net_assigns)
    assigns = net_assigns;
else
    
    error('What kind of file are you trying touse to import assignments? Use a dscalar.nii, pscalar.nii, .txt, or .csv (1 assingment per line.)')
end


%load dconn
disp('Loading dconn and assignments...')
if isnumeric(dconn_cifti_path) ==1
    dconn = dconn_cifti_path;
    clear dconn_cifti_path
else
    dconn_cifti=ciftiopen(dconn_cifti_path,wb_command);
    dconn = single(dconn_cifti.cdata);
end
%subsample dconn



if apply_Zscore_dconn ==1
    newdconn = Zscore_dconn_var(dconn);
else
    newdconn = dconn;
end

clear dconn_cifti dconn

disp('Sorting dconn...')
nonzerounsorted_assigns = assigns~=0;
[sort1,I] = sort(assigns); % get sorted indices;

if exclude_zero_networks ==1
    greys_to_use = sort1~=0;
    networks = unique(nonzeros(assigns));
    adjI = I(greys_to_use);
    I = adjI;
    sorted_dconn1 = newdconn(I,I);
else
    networks = unique(assigns); % get network assingments from template.
    sorted_dconn1 = newdconn(I,I);
end
clear dconn

if plot2dconns ==1
    if isnumeric(dconn_cifti_path2) ==0
    dconn_cifti2=ciftiopen(dconn_cifti_path2,wb_command);
    dconn2 = single(dconn_cifti2.cdata);
    else
      dconn2 = dconn_cifti_path2; 
    end
    clear dconn_cifti2
    if apply_Zscore_dconn ==1
        newdconn2 = Zscore_dconn_var(dconn2);
    else
        newdconn2 = dconn2;
    end
    if use_nets1 ==1
        sorted_dconn2 =  newdconn2(I,I);
    else
        net_assigns2 = ciftiopen(net_assigns2_file,wb_command);
        assigns2 = net_assigns2.cdata;
        [sort2,I] = sort(assigns2); % get sorted indices;
        %networks2 = unique(assigns2); % get network assingments from template.
        disp('Sorting dconn...')
        sorted_dconn2 = newdconn2(I,I);
        clear dconn2
    end
end

f = figure();
set(gcf, 'color','w')

if plot2dconns ==1
    ax1 = subplot(1,3,1);
end
if use_showM ==1
    
else
if downsample_dconn == 1
    if exclude_zero_networks ==1
        assigns_sonly_assigned= assigns(find(assigns~=0));
        assigns_small = assigns_sonly_assigned(1:DS_factor:end);
        
    else
        assigns_small = assigns(1:DS_factor:end);
    end
    
    [sort1small,~] = sort(assigns_small);
    
    sorted_small1 = sorted_dconn1(1:DS_factor:end,1:DS_factor:end);
    imagesc(sorted_small1); hold on;
    insert_net_lines(networks,sort1small,1);
    
else
    imagesc(sorted_dconn1); hold on;
    insert_net_lines(networks,sort1,1)
end
end
set(gca,'FontSize',9)
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.05))
%xlim([0 1]);ylim([0 1]);
%f.PaperPositionMode   = 'auto';
title('Correlation matrix sorted by network','FontSize',9);
if Pos_neg_colormap==1
    load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Positive-Negative_ColorMap.mat','pos_neg_cmap');
    colormap(pos_neg_cmap);
else
colormap jet
end
colorbar;
caxis(caxis_scale);
%f.Position = [100 100 600 600];

if plot2dconns ==1
    ax2=subplot(1,3,2);
    
    if downsample_dconn == 1
        sorted_small2 = sorted_dconn2(1:DS_factor:end,1:DS_factor:end);
        imagesc(sorted_small2); hold on;
        if use_nets1 ==0
            if exclude_zero_networks ==1
                assigns_sonly_assigned= assigns2(find(assigns2~=0));
                assigns_small2 = assigns_sonly_assigned(1:DS_factor:end);
                
            else
                assigns_small2 = assigns2(1:DS_factor:end);
            end
            
            [sort2small,~] = sort(assigns_small2);
            insert_net_lines(networks,sort2small,1);
        else
            
            insert_net_lines(networks,sort1small,1);
        end
        
    else
        imagesc(sorted_dconn2); hold on;
        insert_net_lines(networks,sort2,1)
    end
    
    set(gca,'FontSize',9)
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.05))
    %xlim([0 1]);ylim([0 1]);
    %f.PaperPositionMode   = 'auto';
    title('Correlation matrix sorted by network','FontSize',9);
    colormap jet
    colorbar;
    caxis(caxis_scale)
    %f.Position = [100 100 600 600];
    
    ax3 = subplot(1,3,3);
    if downsample_dconn == 1
        diff_matrix = sorted_small1-sorted_small2;
    else
        diff_matrix = sorted_dconn1-sorted_dconn2;
    end
    imagesc(diff_matrix); hold on;
    insert_net_lines(networks,sort1small,1);
    
    set(gca,'FontSize',9)
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.05))
    title('Matrix Difference','FontSize',9);
    %colormap jet
    load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Positive-Negative_ColorMap.mat','pos_neg_cmap');
    colormap(ax3,pos_neg_cmap);
    colorbar;
    caxis([-0.5 0.5])
    f.Position = [50 100 1300 400];
     disp('saving image...')
    print([image_name '.png'], '-dpng', '-r600')
else
        disp('saving image...')
    print([image_name '.png'], '-dpng', '-r600')
end

if save_non_zeromatrix ==1
    disp('saving...')

    if only_save_abs_diff_dconn ==1
        if exclude_zero_networks ==1
        disp('saving smaller matrix that excludes the zero network assignments (unsorted)...')
        diff_matrix =newdconn(nonzerounsorted_assigns,nonzerounsorted_assigns) - newdconn2(nonzerounsorted_assigns,nonzerounsorted_assigns);
        else
        end
        diff_matrix = abs(diff_matrix);
        save([image_name '_ABS_Diff_dconn1.mat'],'diff_matrix','assigns_sonly_assigned','-v7.3')
    else
        save([image_name '_dconn1.mat'],'sorted_dconn1','assigns_sonly_assigned','-v7.3')
        if plot2dconns ==1
            save([image_name '_dconn1.mat'],'sorted_dconn2','assigns_sonly_assigned','-v7.3')
        end
    end
else
end
disp('Done plotting.')
end

%% insert net lines
function insert_net_lines(networks,sorted_networks,usenet_colors, full_black_lines )
net_start_idx = zeros(1, size(networks,1)); % preallocate for speed
net_end_idx = zeros(1, size(networks,1)); % preallocate for speed
netRGBs = [
    255 0 0;
    0 0 153
    255 255 0
    0 255 0
    13 133 160
    50 50 50
    102 0 204
    102 255 255
    255 128 0
    178 102 153
    0 102 153
    102 255 102
    60 60 251
    200 200 200]/255;
for i=1:size(networks,1)
    netidices = find(sorted_networks == networks(i));
    net_start_idx(i) = min(netidices)-1;
    net_end_idx (i)= max(netidices);
end

    for j = 1:size(networks,1)
        line([net_start_idx(j)-0.5  net_start_idx(j)-0.5], [0.5 size(sorted_networks,1)+0.5],'Color','black','LineWidth',0.5)
        line([0.5 size(sorted_networks,1)+0.5 ],[net_start_idx(j)-0.5  net_start_idx(j)-0.5],'Color','black','LineWidth',0.5)
    end
% 
%     
%     for i = 0:size(H,2)+1
%         line([i-0.5 i-0.5], [0.5 size(H,2)+0.5],'Color','black','LineWidth',2)
%     end
%     for i = 0:size(H,1)+1
%         line([0.5 size(H,1)+0.5],[i-0.5 i-0.5] ,'Color','black','LineWidth',2)
%     end
    

if usenet_colors ==1
    for j = 1:size(networks,1)
        line([net_start_idx(j) net_start_idx(j)], [net_end_idx(j) net_start_idx(j)],'Color',netRGBs(j,:),'LineWidth',2)
        line([net_end_idx(j) net_start_idx(j)], [net_end_idx(j) net_end_idx(j)],'Color',netRGBs(j,:),'LineWidth',2)
        line([net_end_idx(j) net_end_idx(j)], [net_end_idx(j) net_start_idx(j)],'Color',netRGBs(j,:),'LineWidth',2)
        line([net_end_idx(j) net_start_idx(j)], [net_start_idx(j) net_start_idx(j)],'Color',netRGBs(j,:),'LineWidth',2)
    end
else
    for j = 1:size(networks,1)
        line([net_start_idx(j) net_start_idx(j)], [net_end_idx(j) net_start_idx(j)],'Color','black','LineWidth',2)
        line([net_end_idx(j) net_start_idx(j)], [net_end_idx(j) net_end_idx(j)],'Color','black','LineWidth',2)
        line([net_end_idx(j) net_end_idx(j)], [net_end_idx(j) net_start_idx(j)],'Color','black','LineWidth',2)
        line([net_end_idx(j) net_start_idx(j)], [net_start_idx(j) net_start_idx(j)],'Color','black','LineWidth',2)
    end
end
colormap jet
end

%% Zscore transformation
function [newdconn] = Zscore_dconn_var(dconn)
%zscore regions
disp ('sectioning dconn')
LL = dconn(1:29696,1:29696);
LR = dconn(29697:59412,1:29696);
LS = dconn(59413:91282,1:29696);
RL = dconn(1:29696,29697:59412);
RR = dconn(29697:59412,29697:59412);
RS = dconn(59413:91282,29697:59412);
SL = dconn(1:29696,59413:91282);
SR = dconn(29697:59412,59413:91282);
SS = dconn(59413:91282,59413:91282);
disp('resphaping nonants to calcuate Z scores')
ZLL = zscore(reshape(LL,size(LL,1)*size(LL,2),1));
ZLLmat = reshape(ZLL,size(LL,1),size(LL,2)); clear ZLL LL
disp('ZLL');
ZLR = zscore(reshape(LR,size(LR,1)*size(LR,2),1));
ZLRmat = reshape(ZLR,size(LR,1),size(LR,2)); clear ZLR LR
disp('ZLR');
ZLS = zscore(reshape(LS,size(LS,1)*size(LS,2),1));
ZLSmat = reshape(ZLS,size(LS,1),size(LS,2)); clear ZLS LS
disp('ZLS');
ZRL = zscore(reshape(RL,size(RL,1)*size(RL,2),1));
ZRLmat = reshape(ZRL,size(RL,1),size(RL,2)); clear ZRL RL
disp('ZRL');
ZRR = zscore(reshape(RR,size(RR,1)*size(RR,2),1));
ZRRmat = reshape(ZRR,size(RR,1),size(RR,2)); clear ZRR RR
disp('ZRR');
ZRS = zscore(reshape(RS,size(RS,1)*size(RS,2),1));
ZRSmat = reshape(ZRS,size(RS,1),size(RS,2)); clear ZRS RS
disp('ZRS');
ZSL = zscore(reshape(SL,size(SL,1)*size(SL,2),1));
ZSLmat = reshape(ZSL,size(SL,1),size(SL,2)); clear ZSL SL
disp('ZSL');
ZSR = zscore(reshape(SR,size(SR,1)*size(SR,2),1));
ZSRmat = reshape(ZSR,size(SR,1),size(SR,2)); clear ZSR SR
disp('ZSR');
ZSS = zscore(reshape(SS,size(SS,1)*size(SS,2),1));
ZSSmat = reshape(ZSS,size(SS,1),size(SS,2)); clear ZSS SS
disp('ZSS');
disp('rewriting matrix')
newdconn=single(zeros(91282,91282));
newdconn(1:29696,1:29696) = ZLLmat;
newdconn(29697:59412,1:29696) = ZLRmat;
newdconn(59413:91282,1:29696) = ZLSmat;
newdconn(1:29696,29697:59412) = ZRLmat;
newdconn(29697:59412,29697:59412) = ZRRmat;
newdconn(59413:91282,29697:59412) = ZRSmat;
newdconn(1:29696,59413:91282) = ZSLmat;
newdconn(29697:59412,59413:91282) = ZSRmat;
newdconn(59413:91282,59413:91282) = ZSSmat;
clear ZLLmat ZLRmat ZLSmat ZRLmat
clear ZRRmat ZRSmat ZSLmat ZSSmat
clear ZSRmat
end