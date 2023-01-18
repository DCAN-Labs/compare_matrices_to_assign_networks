function plot_surface_mesh
%close all


%% Settings
% outputname_for_cifti_file = 'MSC01_electrode_mappings_BA46';
% load('C:\Users\hermosir\Documents\test_ciftis\MSC01_surfaces\T1w_fsaverage\MSC01_T1w_fsaverage_LR32k_xzy_raw_coordinates.mat','allxyz')
% Lstl='C:\Users\hermosir\Documents\test_ciftis\MSC01_surfaces\T1w_fsaverage\rh.MSC01.L.pial.32k_fs_LR.surf.stl';
% Rstl='C:\Users\hermosir\Documents\test_ciftis\MSC01_surfaces\T1w_fsaverage\rh.MSC01.R.pial.32k_fs_LR.surf.stl';
% path_sulc_dscalar='C:\Users\hermosir\Documents\test_ciftis\MSC01_surfaces\MNINonLinear_fsaverage_LR32K\MSC01.sulc.32k_fs_LR.dscalar.nii';
% network_dscalar='C:\Users\hermosir\Documents\test_ciftis\MSC_to_DCAN\MSC01_half1_to_ADHD315Z_Zscored_recolored.dscalar.nii';

outputname_for_cifti_file = 'PCS0001_electrode_mappings_BA46';
%ABCD files
% load('C:\Users\hermosir\Desktop\TRD_precision_neuromod\sub-PCS0001_T1w_fsaverage_LR32k_xzy_raw_coordinates.mat','allxyz')
% Lstl='C:\Users\hermosir\Desktop\TRD_precision_neuromod\T1w-fsaverage_LR32k-PCS0001.L.pial.32k_fs_LR.surf.stl';
% Rstl='C:\Users\hermosir\Desktop\TRD_precision_neuromod\T1w-fsaverage_LR32k-PCS0001.R.pial.32k_fs_LR.surf.stl';
path_sulc_dscalar='C:\Users\hermosir\Desktop\TRD_precision_neuromod\PCS0001.sulc.32k_fs_LR.dscalar.nii';

%xcp
load('C:\Users\hermosir\Desktop\TRD_precision_neuromod\surface_coords\xcp_pipeline\sub-PCS0001_T1w_fsaverage_LR32k_xzy_raw_coordinates.mat','allxyz');
Lstl='C:\Users\hermosir\Desktop\TRD_precision_neuromod\surface_coords\xcp_pipeline\resampled_sub-PCS0001_ses-01_hemi-L_pial.surf.stl';
Rstl='C:\Users\hermosir\Desktop\TRD_precision_neuromod\surface_coords\xcp_pipeline\resampled_sub-PCS0001_ses-01_hemi-R_pial.surf.stl';

network_dscalar='C:\Users\hermosir\Desktop\TRD_precision_neuromod\TM_networks\ses-01to02_with_ses-01_bad_runs_excluded_smoothed2.55_withoutlier_detection\sub-PCS0001_ses-01to02_task-restMENORDICrmnoisevols_space-fsLR_den-91k_desc-denoisedDilated30mm_bold_SMOOTHED_2.55_bad_runs_excluded_Zscored_recolored.dscalar.nii';



num_stim_paddels = 4;

%ctscan_placements = [25995 27459 54809 56935]; %[ left medial; left lateral; right medial right lateral];
%ctscan_placements = [26100 26885 55673 46458]; rough estimates
ctscan_placements = [26323 27379 57246 55936]; %intial 46 bilateral BA9 placements

ctscan_rotations = [-20 -40 80 -45];
test_vertex(1,1) = ctscan_placements(1);

electrode_array_orient(1,1) = ctscan_rotations(1); %70; % in degrees
test_vertex(2,1) = ctscan_placements(2);
electrode_array_orient(2,1) = ctscan_rotations(2); %180
test_vertex(3,1) = ctscan_placements(3);
electrode_array_orient(3,1) = ctscan_rotations(3); %0
test_vertex(4,1) = ctscan_placements(4);
electrode_array_orient(4,1) = ctscan_rotations(4); %230

%test_vertex = 7998; %try 14533 27078 28155 %39629-Rhemi BA10-28481 BA46:27078:orient 180
%electrode_array_orient = 70; % rotation in degrees

r=14; % paddle lenth is 28mm, so we set the radius of the sphere to be half be half that length.
sulc_depth_threshold =0.25; %typically use 0.35 %larger numbers are more superficial. Recommend to set to 0.15 or greater
planar_radius=14;
plot_l_mesh =1;
plot_r_mesh =1;

electrode_sphere_size = 2.3585; %comes from the length of the diagonla of the contact electrode.

num_electrode_contacts = 8;
electrode_sphere_edge_color = 'none';
electrode_sphere_face_color = [0.7 0.7 0.7];
electrode_sphere_facealpha = 0.2;
run_locally =1;

%% Add dependecies
%add cifti paths
if run_locally ==1
    %Some hardcodes:
    wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
    addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
    addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
    %support_folder='C:\Users\hermosir\Documents\repos\support_folder';
    addpath(genpath('C:\Users\hermosir\Documents\repos\support_folder'));
    load('C:\Users\hermosir\Documents\repos\support_folder\PowerColorMap_wzero.mat')
else
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
    warning('on')
    wb_command=settings.path_wb_c; %path to wb_command
end

netRGBs = [
    255 0 0;
    0 0 153
    255 255 0
    255 255 255
    0 255 0
    255 255 255
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

xyz = allxyz(:,2:4);

figure()
%subplot(1,3,1)
% scatter3(xyz(1:59412,1),xyz(1:59412,2),xyz(1:59412,3),0.5,'k')
% set(gca, 'XLim',[-120 100], 'YLim', [-120 100],'ZLim',[-120 100])
% hold on

if plot_l_mesh ==1
    disp('Importing brain geometric models...')
    model = createpde;
    gm =importGeometry(model,Lstl);
    pdegplot(gm)
    hold on
end
if plot_r_mesh ==1
    model2 = createpde;
    gm2 =importGeometry(model2,Rstl);
    pdegplot(gm2)
    hold on
end
sulc_cii = ft_read_cifti_mod(path_sulc_dscalar);
sulc_depth_extra =sulc_cii.data;

%sulcal depth information do not contain greyordinate to
%vertex maping.
ciimappingscalar = ft_read_cifti_mod(network_dscalar);
mapping_idx = ciimappingscalar.brainstructure ==1 | ciimappingscalar.brainstructure ==2;
TM_networks_full = ciimappingscalar.data;
sulc_depth = sulc_depth_extra(mapping_idx);
TM_networks = TM_networks_full(1:size(sulc_depth,1));
xyz_surface_only = xyz(1:size(sulc_depth),:);


% figure()
% scatter3(xyz_surface_only(:,1),xyz_surface_only(:,2),xyz_surface_only(:,3),0.5,'k')
% set(gca, 'XLim',[-120 100], 'YLim', [-120 100],'ZLim',[-120 100])
sulcal_indx = sulc_depth > sulc_depth_threshold;

%unique_nets = unique(TM_networks);
netcolors_persulc = zeros(size(sulc_depth,1),3);

%for i = 1:size(unique_nets,1)
for i = 1:size(TM_networks,1)
    this_nets_value = TM_networks(i);
    netcolors_persulc(i,:) = netRGBs(this_nets_value,:);
    %                     this_nets_idxs = TM_networks == i;
    %                     netcolors_persulc(this_nets_idxs,:) =netRGBs(unique_nets(i),:);
end
net_colors_surf = netcolors_persulc(sulcal_indx,1:3);
scatter3(xyz_surface_only(sulcal_indx,1),xyz_surface_only(sulcal_indx,2),xyz_surface_only(sulcal_indx,3),8, net_colors_surf,'filled')
set(gca,'Color','k')

%% plot sphere

for elec = 1:num_stim_paddels
    %r=20;
    %define a 20x20 sphere.(20 segments from pole to pole and
    %20 segments around the equator).
    [X,Y,Z] = sphere(20);
    [sphere_all]=xyz_surface_only(test_vertex(elec,1),:);
    X0= sphere_all(1);
    Y0= sphere_all(2);
    Z0=sphere_all(3);
    
    x=X*r + X0;
    y=Y*r + Y0;
    z=Z*r + Z0;
    lightgrey=0.2*[1 1 1];
    
    %Plot sphere where the the diameter of the length of the array  
    surface(x,y,z,'FaceColor','none','EdgeColor',lightgrey)
    
    hold on
    
    gyri_points = xyz_surface_only(sulcal_indx,1:3);
    gyri_points_orig_mappings{elec,1} = find(sulcal_indx ==1);
    
    disp('calculating between vertices distances...')
    for i = 1:size(gyri_points,1)
        distance_vec(i,1) = sqrt((sphere_all(1) - gyri_points(i,1))^2 + (sphere_all(2) - gyri_points(i,2))^2 + (sphere_all(3) - gyri_points(i,3))^2 );
    end
    all_distance_vec{1,elec} = distance_vec;
end

figure()
%subplot(1,3,2)
grey_surfacepoints = ones(size(xyz_surface_only,1),3)*0.3;
scatter3(xyz_surface_only(sulcal_indx,1),xyz_surface_only(sulcal_indx,2),xyz_surface_only(sulcal_indx,3),6, grey_surfacepoints(sulcal_indx,1:3),'filled')
%Background color of target plot
set(gca,'Color','w'); hold on

for elec = 1:num_stim_paddels
    distance_vec = cell2mat(all_distance_vec(1,elec));
    all_within_sphere_idx{elec,1} = find(distance_vec<r); %find the grayordinates that are within the sphere radius.
    within_sphere_idx = find(distance_vec<r);
    %net_colors_surf(within_sphere_idx);
    
    %Plot grayordinates within radius
    scatter3(gyri_points(within_sphere_idx,1),gyri_points(within_sphere_idx,2),gyri_points(within_sphere_idx,3),12, net_colors_surf(within_sphere_idx,1:3),'filled')
    %Plot center point
    scatter3(xyz_surface_only(test_vertex(elec,1),1),xyz_surface_only(test_vertex(elec,1),2),xyz_surface_only(test_vertex(elec,1),3),50, netcolors_persulc(test_vertex(elec,1),1:3),'filled')
    set(gca,'Color','k')
    hold on;
    
    %within_sphere_xyz(elec,1) = [gyri_points(within_sphere_idx(elec,1),1) gyri_points(within_sphere_idx(elec,1),2) gyri_points(within_sphere_idx(elec,1),3)];
    
    [planar_idxs] = find(distance_vec<planar_radius); %find the grayordinates that are within the sphere radius.
    planar_xyzs{elec} = [gyri_points(planar_idxs,1) gyri_points(planar_idxs,2) gyri_points(planar_idxs,3)];
end

axis equal
box on


%% Get normal plane
%get normals
for elec = 1:num_stim_paddels
    
    [n_2,V_2,p_2] = affine_fit(cell2mat(planar_xyzs(elec)));
    n_2_all{elec} = n_2;
    V_2_all{elec} = V_2;
    p_2_all{elec} = p_2;
    
    n_2 = cell2mat(n_2_all(elec));
    V_2 = cell2mat(V_2_all(elec));
    p_2 = cell2mat(p_2_all(elec));
    
    n_2_native = [-sqrt(2)/2 0 sqrt(2)/2];
    %quiver3(p_2(1),p_2(2),p_2(3),n_2(1)/3,n_2(2)/3,n_2(3)/3,'r','linewidth',10,'Autoscale',1)
    quiver3(p_2(1),p_2(2),p_2(3),n_2(1)*5,n_2(2)*5,n_2(3)*5,'r','linewidth',5,'Autoscale',1)
    [S1,S2] = meshgrid([-10 0 10]);
    %generate the point coordinates
    X = p_2(1)+[S1(:) S2(:)]*V_2(1,:)';
    Y = p_2(2)+[S1(:) S2(:)]*V_2(2,:)';
    Z = p_2(3)+[S1(:) S2(:)]*V_2(3,:)';
    
    %plot the plane
    %surf(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),'facecolor','red','facealpha',0.3);
    
    %%build circular mesh
    M = 8 ; %Number of rings
    N = 36 ; %number of steps
    R1 = 0 ; % inner radius
    R2 = r ;  % outer radius
    nR = linspace(R1,R2,M) ;
    nT = linspace(0,2*pi,N) ;
    %nT = pi/180*(0:NT:theta) ;
    %[meshR, T] = meshgrid(nR,nT) ;
    [meshR, T] = meshgrid(linspace(0,2*pi,N),linspace(R1,R2,M));
    
    % Convert grid to cartesian coordintes
    Xvec1=T.*cos(meshR);
    Xvec2=Xvec1(:); %put the cartsian transformed data into a vector.
    Yvec1=T.*sin(meshR);
    Yvec2=Yvec1(:);
    %Z=Xvec1.*2;
    % Xvec1=cos(T); Xvec2=Xvec1(:); %put the cartsian transformed data into a vector.
    % Yvec1=sin(T); Yvec2=Yvec1(:);
    
    X = p_2(1)+ [Xvec2 Yvec2]*V_2(1,:)';
    Y = p_2(2)+ [Xvec2 Yvec2]*V_2(2,:)';
    Z = p_2(3)+ zeros(size(Xvec2,1),2)*V_2(3,:)'; % The X-Y plane is linear, you could deform it here based on the radius of the brain.
    Z = p_2(3)+ [Xvec2 Yvec2]*V_2(3,:)';
    %[m,n]=size(X);
    %surf(X*r,Y*r,Z*r,'FaceAlpha',0.5,'EdgeColor','w')
    %surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','w')
    X_reshaped = reshape(X,M,N);
    Y_reshaped = reshape(Y,M,N);
    Z_reshaped = reshape(Z,M,N);
    surf(X_reshaped,Y_reshaped,Z_reshaped, 'FaceAlpha',0,'EdgeColor',[0.7 0.7 0.7])
end

xlabel('x');
ylabel('y');
zlabel('z');
axis equal

%% try to plot a cube
% edges = [7,r*2,2];  %assume array with is 28 * 7 * 2
% %origin = [0,0,0];
% %[edges,origin,alpha,clr,rotate_vec,origin_point] = deal(inArgs{:});
% clr = [1,0,0];
% alpha = 0.3;
% origin = [0 0 0];
% origin(1) = origin(1)-(edges(1)/2);
% origin(2) = origin(2)-(edges(2)/2);
% origin(3) = origin(3)-(edges(3)/2);

% XYZ = { ...
%   [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
%   [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
%   [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
%   [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
%   [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
%   [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
%   };
%
% XYZ = mat2cell(...
%   cellfun( @(x,y,z) x*y+z , ...
%     XYZ , ...
%     repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
%     repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
%     'UniformOutput',false), ...
%   6,[1 1 1]);
%
% % The origin offset has already been done, so the only the linear
% % tranformation remains.
% XYZ_prime = XYZ;
%     for k =1:size(XYZ{1},1)
%             Xpts = XYZ{1,1}{k,1}';
%             Ypts = XYZ{1,2}{k,1}';
%             Zpts = XYZ{1,3}{k,1}';
%             %Xprime = [Xpts Ypts]*V_2(1,:)' -(edges(1)/2) + p_2(1);
%             %Yprime = [Xpts Ypts]*V_2(2,:)' -(edges(2)/2) + p_2(2);
%             %Zprime = [Xpts Zpts]*V_2(3,:)' -(edges(3)/2) + p_2(3);
%             Xprime = [Xpts Ypts]*V_2(1,:)' + p_2(1);
%             Yprime = [Xpts Ypts]*V_2(2,:)' + p_2(2);
%             Zprime = [Xpts Ypts]*V_2(3,:)' + p_2(3);
%             XYZ_prime{1,1}{k,1} = Xprime';
%             XYZ_prime{1,2}{k,1} = Yprime';
%             XYZ_prime{1,3}{k,1} = Zprime';
%     end
%plot electrode orientation
%rg = cellfun(@patch,XYZ_prime{1},XYZ_prime{2},XYZ_prime{3},repmat({clr},6,1),repmat({'FaceAlpha'},6,1),repmat({alpha},6,1));

%figure()
%leed array spacing
all_electrode_origins = cell(num_stim_paddels,1);

for elec = 1:num_stim_paddels
    leed1 = polyshape([0 2.5 2.5 0],[0 0 4 4]);
    leed2 = polyshape([4.5 7 7 4.5],[3 3 7 7]);
    leed3 = polyshape([0 2.5 2.5 0],[7 7 11 11]);
    leed4 = polyshape([4.5 7 7 4.5],[10 10 14 14]);
    leed5 = polyshape([0 2.5 2.5 0],[14 14 18 18]);
    leed6 = polyshape([4.5 7 7 4.5],[17 17 21 21]);
    leed7 = polyshape([0 2.5 2.5 0],[21 21 25 25]);
    leed8 = polyshape([4.5 7 7 4.5],[24 24 28 28]);
    
    array_center = [3.5 14]; % if the array started with the corner at (0,0) this would be the center, however, the center needs to moved to (0,0);
    lead_origin = [0,0];
    
    leed1_origin_adjusted_shape = polyshape([0 2.5 2.5 0]-array_center(1),[0 0 4 4]-array_center(2));
    leed2_origin_adjusted_shape = polyshape([4.5 7 7 4.5]-array_center(1),[3 3 7 7]-array_center(2));
    leed3_origin_adjusted_shape = polyshape([0 2.5 2.5 0]-array_center(1),[7 7 11 11]-array_center(2));
    leed4_origin_adjusted_shape = polyshape([4.5 7 7 4.5]-array_center(1),[10 10 14 14]-array_center(2));
    leed5_origin_adjusted_shape = polyshape([0 2.5 2.5 0]-array_center(1),[14 14 18 18]-array_center(2));
    leed6_origin_adjusted_shape = polyshape([4.5 7 7 4.5]-array_center(1),[17 17 21 21]-array_center(2));
    leed7_origin_adjusted_shape = polyshape([0 2.5 2.5 0]-array_center(1),[21 21 25 25]-array_center(2));
    leed8_origin_adjusted_shape = polyshape([4.5 7 7 4.5]-array_center(1),[24 24 28 28]-array_center(2));
    
    leed1prime = rotate(leed1_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed2prime = rotate(leed2_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed3prime = rotate(leed3_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed4prime = rotate(leed4_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed5prime = rotate(leed5_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed6prime = rotate(leed6_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed7prime = rotate(leed7_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    leed8prime = rotate(leed8_origin_adjusted_shape,electrode_array_orient(elec,1),lead_origin);
    
    leed1prime_verts = leed1prime.Vertices;
    X1 = p_2_all{elec}(1)+ leed1prime_verts*V_2_all{elec}(1,:)';
    Y1 = p_2_all{elec}(2)+ leed1prime_verts*V_2_all{elec}(2,:)';
    Z1 = p_2_all{elec}(3)+ leed1prime_verts*V_2_all{elec}(3,:)';
    % get mean for electrode sphere
    X1_rotated_origin = mean(X1); all_electrode_origins{elec}(1,1) = X1_rotated_origin;
    Y1_rotated_origin = mean(Y1); all_electrode_origins{elec}(1,2) = Y1_rotated_origin;
    Z1_rotated_origin = mean(Z1); all_electrode_origins{elec}(1,3) = Z1_rotated_origin;
    
    leed2prime_verts = leed2prime.Vertices;
    X2 = p_2_all{elec}(1)+ leed2prime_verts*V_2_all{elec}(1,:)';
    Y2 = p_2_all{elec}(2)+ leed2prime_verts*V_2_all{elec}(2,:)';
    Z2 = p_2_all{elec}(3)+ leed2prime_verts*V_2_all{elec}(3,:)';
    X2_rotated_origin = mean(X2); all_electrode_origins{elec}(2,1) = X2_rotated_origin;
    Y2_rotated_origin = mean(Y2); all_electrode_origins{elec}(2,2) = Y2_rotated_origin;
    Z2_rotated_origin = mean(Z2); all_electrode_origins{elec}(2,3) = Z2_rotated_origin;
    
    leed3prime_verts = leed3prime.Vertices;
    X3 = p_2_all{elec}(1)+ leed3prime_verts*V_2_all{elec}(1,:)';
    Y3 = p_2_all{elec}(2)+ leed3prime_verts*V_2_all{elec}(2,:)';
    Z3 = p_2_all{elec}(3)+ leed3prime_verts*V_2_all{elec}(3,:)';
    X3_rotated_origin = mean(X3); all_electrode_origins{elec}(3,1) = X3_rotated_origin;
    Y3_rotated_origin = mean(Y3); all_electrode_origins{elec}(3,2) = Y3_rotated_origin;
    Z3_rotated_origin = mean(Z3); all_electrode_origins{elec}(3,3) = Z3_rotated_origin;
    
    leed4prime_verts = leed4prime.Vertices;
    X4 = p_2_all{elec}(1)+ leed4prime_verts*V_2_all{elec}(1,:)';
    Y4 = p_2_all{elec}(2)+ leed4prime_verts*V_2_all{elec}(2,:)';
    Z4 = p_2_all{elec}(3)+ leed4prime_verts*V_2_all{elec}(3,:)';
    X4_rotated_origin = mean(X4); all_electrode_origins{elec}(4,1) = X4_rotated_origin;
    Y4_rotated_origin = mean(Y4); all_electrode_origins{elec}(4,2) = Y4_rotated_origin;
    Z4_rotated_origin = mean(Z4); all_electrode_origins{elec}(4,3) = Z4_rotated_origin;
    
    leed5prime_verts = leed5prime.Vertices;
    X5 = p_2_all{elec}(1)+ leed5prime_verts*V_2_all{elec}(1,:)';
    Y5 = p_2_all{elec}(2)+ leed5prime_verts*V_2_all{elec}(2,:)';
    Z5 = p_2_all{elec}(3)+ leed5prime_verts*V_2_all{elec}(3,:)';
    X5_rotated_origin = mean(X5); all_electrode_origins{elec}(5,1) = X5_rotated_origin;
    Y5_rotated_origin = mean(Y5); all_electrode_origins{elec}(5,2) = Y5_rotated_origin;
    Z5_rotated_origin = mean(Z5); all_electrode_origins{elec}(5,3) = Z5_rotated_origin;
    
    leed6prime_verts = leed6prime.Vertices;
    X6 = p_2_all{elec}(1)+ leed6prime_verts*V_2_all{elec}(1,:)';
    Y6 = p_2_all{elec}(2)+ leed6prime_verts*V_2_all{elec}(2,:)';
    Z6 = p_2_all{elec}(3)+ leed6prime_verts*V_2_all{elec}(3,:)';
    X6_rotated_origin = mean(X6); all_electrode_origins{elec}(6,1) = X6_rotated_origin;
    Y6_rotated_origin = mean(Y6); all_electrode_origins{elec}(6,2) = Y6_rotated_origin;
    Z6_rotated_origin = mean(Z6); all_electrode_origins{elec}(6,3) = Z6_rotated_origin;
    
    leed7prime_verts = leed7prime.Vertices;
    X7 = p_2_all{elec}(1)+ leed7prime_verts*V_2_all{elec}(1,:)';
    Y7 = p_2_all{elec}(2)+ leed7prime_verts*V_2_all{elec}(2,:)';
    Z7 = p_2_all{elec}(3)+ leed7prime_verts*V_2_all{elec}(3,:)';
    X7_rotated_origin = mean(X7); all_electrode_origins{elec}(7,1) = X7_rotated_origin;
    Y7_rotated_origin = mean(Y7); all_electrode_origins{elec}(7,2) = Y7_rotated_origin;
    Z7_rotated_origin = mean(Z7); all_electrode_origins{elec}(7,3) = Z7_rotated_origin;
    
    leed8prime_verts = leed8prime.Vertices;
    X8 = p_2_all{elec}(1)+ leed8prime_verts*V_2_all{elec}(1,:)';
    Y8 = p_2_all{elec}(2)+ leed8prime_verts*V_2_all{elec}(2,:)';
    Z8 = p_2_all{elec}(3)+ leed8prime_verts*V_2_all{elec}(3,:)';
    X8_rotated_origin = mean(X8); all_electrode_origins{elec}(8,1) = X8_rotated_origin;
    Y8_rotated_origin = mean(Y8); all_electrode_origins{elec}(8,2) = Y8_rotated_origin;
    Z8_rotated_origin = mean(Z8); all_electrode_origins{elec}(8,3) = Z8_rotated_origin;
    
    fill3(X1,Y1,Z1,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X2,Y2,Z2,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X3,Y3,Z3,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X4,Y4,Z4,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X5,Y5,Z5,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X6,Y6,Z6,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X7,Y7,Z7,[0.9 0.9 0.9],'FaceAlpha',0.5)
    fill3(X8,Y8,Z8,[0.9 0.9 0.9],'FaceAlpha',0.5)
    %A =plotcube(X,Y,Z);
    
    scatter3(X1_rotated_origin,Y1_rotated_origin,Z1_rotated_origin,4,'w','filled')
    scatter3(X2_rotated_origin,Y2_rotated_origin,Z2_rotated_origin,4,'w','filled')
    scatter3(X3_rotated_origin,Y3_rotated_origin,Z3_rotated_origin,4,'w','filled')
    scatter3(X4_rotated_origin,Y4_rotated_origin,Z4_rotated_origin,4,'w','filled')
    scatter3(X5_rotated_origin,Y5_rotated_origin,Z5_rotated_origin,4,'w','filled')
    scatter3(X6_rotated_origin,Y6_rotated_origin,Z6_rotated_origin,4,'w','filled')
    scatter3(X7_rotated_origin,Y7_rotated_origin,Z7_rotated_origin,4,'w','filled')
    scatter3(X8_rotated_origin,Y8_rotated_origin,Z8_rotated_origin,4,'w','filled')
    
    %create spheres at electrodes
    [sphere1X,sphere1Y,sphere1Z] = sphere(15);
    sph1x=sphere1X*electrode_sphere_size + X1_rotated_origin;
    sph1y=sphere1Y*electrode_sphere_size + Y1_rotated_origin;
    sph1z=sphere1Z*electrode_sphere_size + Z1_rotated_origin;
    surface(sph1x,sph1y,sph1z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph2x=sphere1X*electrode_sphere_size + X2_rotated_origin;
    sph2y=sphere1Y*electrode_sphere_size + Y2_rotated_origin;
    sph2z=sphere1Z*electrode_sphere_size + Z2_rotated_origin;
    surface(sph2x,sph2y,sph2z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph3x=sphere1X*electrode_sphere_size + X3_rotated_origin;
    sph3y=sphere1Y*electrode_sphere_size + Y3_rotated_origin;
    sph3z=sphere1Z*electrode_sphere_size + Z3_rotated_origin;
    surface(sph3x,sph3y,sph3z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph4x=sphere1X*electrode_sphere_size + X4_rotated_origin;
    sph4y=sphere1Y*electrode_sphere_size + Y4_rotated_origin;
    sph4z=sphere1Z*electrode_sphere_size + Z4_rotated_origin;
    surface(sph4x,sph4y,sph4z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph5x=sphere1X*electrode_sphere_size + X5_rotated_origin;
    sph5y=sphere1Y*electrode_sphere_size + Y5_rotated_origin;
    sph5z=sphere1Z*electrode_sphere_size + Z5_rotated_origin;
    surface(sph5x,sph5y,sph5z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph6x=sphere1X*electrode_sphere_size + X6_rotated_origin;
    sph6y=sphere1Y*electrode_sphere_size + Y6_rotated_origin;
    sph6z=sphere1Z*electrode_sphere_size + Z6_rotated_origin;
    surface(sph6x,sph6y,sph6z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    sph7x=sphere1X*electrode_sphere_size + X7_rotated_origin;
    sph7y=sphere1Y*electrode_sphere_size + Y7_rotated_origin;
    sph7z=sphere1Z*electrode_sphere_size + Z7_rotated_origin;
    surface(sph7x,sph7y,sph7z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
    [sphere1X,sphere1Y,sphere1Z] = sphere(15);
    sph8x=sphere1X*electrode_sphere_size + X8_rotated_origin;
    sph8y=sphere1Y*electrode_sphere_size + Y8_rotated_origin;
    sph8z=sphere1Z*electrode_sphere_size + Z8_rotated_origin;
    surface(sph8x,sph8y,sph8z,'FaceColor',electrode_sphere_face_color,'EdgeColor',electrode_sphere_edge_color,'FaceAlpha',electrode_sphere_facealpha)
    hold on
    
end

for elec = 1:num_stim_paddels
    for k = 1:num_electrode_contacts
        for i = 1:size(gyri_points,1)
            distance_vec(i,k) = sqrt((all_electrode_origins{elec}(k,1) - gyri_points(i,1))^2 + (all_electrode_origins{elec}(k,2) - gyri_points(i,2))^2 + (all_electrode_origins{elec}(k,3) - gyri_points(i,3))^2 );
        end
        distance_vec_per_electrode{elec,1}= distance_vec;
    end
end

clear within_sphere_idx %used earlier
for elec = 1:num_stim_paddels
    for k = 1:num_electrode_contacts
        %within_sphere_idx = find(distance_vec(:,k)<electrode_sphere_size);
        within_sphere_idx{elec,1} = find(distance_vec_per_electrode{elec,1}(:,k)<electrode_sphere_size);
        within_sphere_gyralindices{elec,1}{k,1} = within_sphere_idx{elec,1};
    end
end

sulcal_TM_networks = TM_networks(sulcal_indx);

for elec = 1:num_stim_paddels
    for k = 1:num_electrode_contacts
        areinsphere = within_sphere_gyralindices{elec,1}{k,1};
        electrode_nets{elec,1}{k,1} = sulcal_TM_networks(areinsphere);
        electrode_nets_mode{elec,1}(k,1) = mode(electrode_nets{elec,1}{k,1});
        electrode_nets_RGBcolors{elec,1}{k,1} = net_colors_surf(areinsphere,:);
        this_net_mode_RGBcolors{elec,1}(k,1) = mode(net_colors_surf(areinsphere));
    end
    
end

% plot(leed1prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed2prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed3prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed4prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed5prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed6prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed7prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
% plot(leed8prime,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.8); hold on;
f  = figure();
for elec = 1:num_stim_paddels
    
    if mod(elec,2)==0 % iseven
        tranlation_vec = [12,(elec*20)-40];
    else %elec is odd
        if elec ==1
            tranlation_vec = [0,0];
        else
            tranlation_vec = [0,((elec+1)*20)-40];
        end
    end
    
    this_elec_lead1 = translate(leed1,tranlation_vec);
    this_elec_lead2 = translate(leed2,tranlation_vec);
    this_elec_lead3 = translate(leed3,tranlation_vec);
    this_elec_lead4 = translate(leed4,tranlation_vec);
    this_elec_lead5 = translate(leed5,tranlation_vec);
    this_elec_lead6 = translate(leed6,tranlation_vec);
    this_elec_lead7 = translate(leed7,tranlation_vec);
    this_elec_lead8 = translate(leed8,tranlation_vec);
    
    
    try
        plot(this_elec_lead1,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(1,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead1,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead2,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(2,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead2,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead3,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(3,1),:),'FaceAlpha',1); hold on;
        
    catch
        plot(this_elec_lead3,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead4,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(4,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead4,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead5,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(5,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead5,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead6,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(6,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead6,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead7,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(7,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead7,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    try
        plot(this_elec_lead8,'FaceColor',netRGBs(electrode_nets_mode{elec,1}(8,1),:),'FaceAlpha',1); hold on;
    catch
        plot(this_elec_lead8,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1); hold on;
    end
    
    if mod(elec,2)==0
        text_location_x = ((this_elec_lead2.Vertices(4,1) - this_elec_lead1.Vertices(1,1))/2) + 7;
    else
        text_location_x = ((this_elec_lead2.Vertices(4,1) - this_elec_lead1.Vertices(1,1))/2) -4;
    end
    text_location_y = this_elec_lead1.Vertices(1,2)-5;
    
    text(text_location_x,text_location_y,['Array ' num2str(elec)])
    
end


axis equal
xlim([-5 (max(tranlation_vec(1)+12))]);
ylim([-10 (max(tranlation_vec(2)+30))]);
set(gcf,'Position',[1000 100 150 450]);
set(gca,'visible','off');

 %% remap the networks to the orginal grayordinate indices.
for elec = 1:num_stim_paddels
    export_dscalar = zeros(size(TM_networks_full,1),1);
    gyral_networks = TM_networks_full(gyri_points_orig_mappings{elec,1});
    
    for m = 1:num_electrode_contacts
        this_electrodes_sulc_indxs = within_sphere_gyralindices{elec,1}{m};
        this_electrodes_greymappings = gyri_points_orig_mappings{elec,1}(this_electrodes_sulc_indxs);
        export_dscalar(this_electrodes_greymappings) = gyral_networks(this_electrodes_sulc_indxs);
    end
    
    saving_Cii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\MSC_to_DCAN\MSC01_half1_to_ADHD315Z_Zscored_recolored.dscalar.nii',wb_command);
    saving_Cii.cdata = export_dscalar;
    ciftisave(saving_Cii,[outputname_for_cifti_file '_array' num2str(elec) '.dscalar.nii'],wb_command);
end


% r = cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
%   repmat({clr},6,1),...
%   repmat({'FaceAlpha'},6,1),...
%   repmat({alpha},6,1)...
%   );

% model = createpde;
% gm =importGeometry(model,'C:\Users\hermosir\Downloads\lamitrode44C_v1.stl');
% [F] = stlread('C:\Users\hermosir\Downloads\lamitrode44C_v1.stl');
% f = pdegplot(gm,'FaceLabels','on');
% patch(F.Points(:,1),F.Points(:,2),F.Points(:,3),'FaceColor',[0.8 0.8 1])
% translate(g,10);

%
% figure()
% drawarrayplot(modes) %modes should be 2 by 4 array with each network assignment in each array.

% normals = pcnormals(within_sphere_xyz);
% p = mean(within_sphere_xyz,1);

% %The samples are reduced:
% R = bsxfun(@minus,within_sphere_xyz,p);
% %Computation of the principal directions if the samples cloud
% [V,D] = eig(R'*R);
% [d,idx_sorted] = sort(diag(D));
% % sort eigenvalues matrix
% D_sorted = D(idx_sorted,idx_sorted);
% % sort eigenvectors
% V_sorted = V(:,idx_sorted);
% %Extract the output from the eigenvectors
% n = V(:,1);
% V = V(:,2:end);


disp('Done.')

end