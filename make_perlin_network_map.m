function [output_dscalar_name] = make_perlin_network_map

close all
n = 500;
m = 500;
reset(RandStream.getGlobalStream,sum(100*clock)); % added to ensure randomness of frame sampling during parellization -RH 02/22/2024
v=num2str(randi([1 10000000]));

%n = 64; %use a smaller matrix to debug faster
%m = 64;
num_nets = 16;
%load('C:\Users\hermosir\Documents\repos\support_folder\PowerColorMap_wzero.mat');
load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/PowerColorMap_wzero.mat');
load('/home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/PowerColorMap.mat')
%load('C:\Users\hermosir\Documents\repos\support_folder\PowerColorMap.mat');
load('/home/rando149/shared/projects/Polyvertexscore/parcel_probability_map.mat');
%load the parcel file to get the RGB colors;
load('/home/faird/shared/code/internal/utilities/distance-matrix/Conte69_flatmap_and_vol_91282_xyz_coord.mat','allxyz')

% provide a dscalar to make a random map from
dscalar_file='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii';

wb_command = '/home/faird/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command';
% addpath(genpath('/home/faird/shared/code/external/utilities/gifti/gifti-1.6'));
% addpath(genpath('/home/faird/shared/code/external/utilities/Matlab_CIFTI'));
addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/'));
addpath(genpath('/home/faird/shared/code/internal/utilities/Matlab_CIFTI/'));
addpath(genpath('/home/miran045/shared/code/internal/utilities/CIFTI/'));
addpath(genpath('/home/miran045/shared/code/internal/utilities/gifti'));


% for i = 1:size(parcel_file.parcel,2)
%     net_colors((parcel_file.parcel(i).ix),:) =  repmat(parcel_file.parcel(i).RGB,size(parcel_file.parcel(i).ix,1),1);
% end
ax =figure; %1

for j=1:num_nets
    im = zeros(n, m);
    all_im(j,:,:) = perlin_noise(im);
    
    subplot(4,5,j)
    imagesc(squeeze(all_im(j,:,:))); colormap jet;
   
    %use x and y to offset the imagesc axis.
    x=[n/2*-1 n/2];
    y=[m/2*-1 m/2];
        pics_coords_x=x(1):x(end)-1;
    pics_coords_y=y(1):y(end)-1;
    
    for xx=1:(x(end)-1-x(1))
        for yy=1:(y(end)-1-y(1))
            pix_coords{xx,yy}=[pics_coords_x(xx), pics_coords_y(yy)];
        end
    end
    imagesc(x,y,squeeze(all_im(j,:,:))); colormap jet;
    set(gca,'YDir','normal')
end
colormap(jet);
% values are expressed as positive and negative values. so, we have to
% raise all the values up to make them all positive.

all_im =all_im+abs((min(all_im(:))))+1;

all_im(4,:,:) = zeros(n, m);
all_im(6,:,:) = zeros(n, m);


[~,B]=max(all_im);
B= squeeze(B);
B_flipped=flip(B,2);

%subplot(3,5,15) %plot on to the last spot.
figure() %2
imagesc(x,y,B);
set(gca,'YDir','normal')
%load('C:\Users\hermosir\Documents\repos\support_folder\PowerColorMap.mat');
colormap(mymap)
%colormap(ax(15),jet)

figure()
%load('Conte69_flatmap_and_vol_91282_xyz_coord.mat')
%scatter(allxyz(1:59412,2),allxyz(1:59412,3));
imagesc(x,y,B);
set(gca,'YDir','normal')
hold on
colormap(mymap)
scatter(allxyz(1:29696,2),allxyz(1:29696,3),'k','.');
set(gca,'Color','k')
L_mymesh=[allxyz(1:29696,2),allxyz(1:29696,3)];
R_mymesh=[allxyz(29697:59412,2),allxyz(29697:59412,3)];


[X,Y] = meshgrid(pics_coords_x,pics_coords_y);
% NOTE: Be careful with this reshape. The index referencing of B, assumes a
% specific order of the indexing.  best to use square dimensions.
X_vec=reshape(X,(size(X,1)*(size(X,2))),1); 
Y_vec=reshape(Y,(size(Y,1)*(size(Y,2))),1);

cart_pairs= [X_vec,Y_vec];
[DL,IL]=pdist2(cart_pairs,L_mymesh,'euclidean','Smallest',1);
[DR,IR]=pdist2(cart_pairs,R_mymesh,'euclidean','Smallest',1);

assignsL =B(IL); 
assignsR =B_flipped(IR);
all_rgbs_L=zeros(size(assignsL,2),4);
all_rgbs_R=zeros(size(assignsR,2),4);

for m=1:size(assignsL,2)
    thisgrays_parcel_row = find([parcel.power_val]==assignsL(m));
    all_rgbs_L(m,:) = parcel(thisgrays_parcel_row).RGB;
end
all_rgbs_L(:,4) = [];

figure();
%scatter(allxyz(1:29696,2),allxyz(1:29696,3),[],all_rgbs_L,'.');
scatter(allxyz(1:29696,2),allxyz(1:29696,3),[],assignsL,'.'); % easier way to color the scatterplot
set(gca,'Color',[0.1 0.1 0.1])
colormap(mymap)
        xlims=[x(1) x(end)-1];
        ylims=[y(1) y(end)-1];
        xlim(xlims);
        ylim(ylims);

for m=1:size(assignsR,2)
    thisgrays_parcel_row = find([parcel.power_val]==assignsR(m));
    all_rgbs_R(m,:) = parcel(thisgrays_parcel_row).RGB;
end
all_rgbs_R(:,4) = [];
figure();
%scatter(allxyz(29697:59412,2),allxyz(29697:59412,3),[],all_rgbs_R,'.');
scatter(allxyz(29697:59412,2),allxyz(29697:59412,3),[],assignsR,'.'); % easier way to color the scatterplot
colormap(mymap)
set(gca,'Color',[0.1 0.1 0.1])
        xlims=[x(1) x(end)-1];
        ylims=[y(1) y(end)-1];
        xlim(xlims);
        ylim(ylims);        
        
% for i=1:numNodes
%     for j=1:size(pix_coords_x)*size(pix_coords_y);
%   distances(i) = sqrt((sphere_all(1) - gyri_points(i,1))^2 + (sphere_all(2) - gyri_points(i,2))^2 + (sphere_all(3) - gyri_points(i,3))^2 );
%         
%     end
% end
cii = ciftiopen(dscalar_file,wb_command);
cii.cdata(1:29696)=assignsL;
cii.cdata(29697:59412)=assignsR;
ciftisave(cii,['Random_map' num2str(v) '.dscalar.nii'],wb_command);
%clean the file
output_dscalar_name = clean_dscalars_by_size(['Random_map' num2str(v) '.dscalar.nii'],[],[],[],[],30,[],0,0,1,0);
disp('Done.')

end


function im = perlin_noise(im)

    [n, m] = size(im);
    i = 0;
    w = sqrt(n*m);
    disp('Converging noise...')
    
    while w > 3
        tic
        %if i>=2 
        i = i + 1;
        %else
        %i = i + 2;            
        %end
        d = interp2(randn(n, m), i-1, 'spline');
        im = im + i * d(1:n, 1:m);
        w = w - ceil(w/2 - 1);
        toc
    end
end