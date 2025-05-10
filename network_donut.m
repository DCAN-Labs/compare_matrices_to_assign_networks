function network_donut(numdat, sorted_labels_a, sorted_labels_b, chart_title,sort_mat,Eidx, Fidx, isanet_in_current_setC, isanet_in_current_setD,skip_plotting,add_label_percentages,donut_label_font_size)

%Some hardcodes:
% wb_command = ('C:\Users\hermosir\Desktop\workbench\bin_windows64\wb_command');
% addpath(genpath('C:\Users\hermosir\Documents\repos\HCP_MATLAB'));
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\utilities')
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\gifti')
% addpath('C:\Users\hermosir\Documents\repos\MSCcodebase-master\Utilities\read_write_cifti\fileio')
%add_label_percentages = 1;
label_threshold = 5; %only label networks greater than 1%  setting this to zero will label all networks (and probably look messy).
%OPTIONS
% close all
% check_surface =1;
% check_subcortical =1;
% sort_matrix =0;
% plot_subcort =1;
%
% all_labels = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SML','AUD', 'Tpole', 'MTL','PMN','PON'};
% possible_net_nums = [1 2 3  5  7 8 9 10 11 12 13 14 15 16];

% %% Load data
% Ccii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\Unique_5minute_Control_mode_paths.dscalar.nii',wb_command);
% C = Ccii.cdata;
% if check_subcortical ==1
%     C_info = ft_read_cifti_mod('C:\Users\hermosir\Documents\test_ciftis\Unique_5minute_Control_mode_paths.dscalar.nii');
%     greybrainstructs_C = C_info.brainstructure(find(C_info.brainstructure >0));
%     surfacegreys_C = size(find(C_info.brainstructure <3 & C_info.brainstructure >0));
%     C_surf = C(1:surfacegreys_C,1);
%     C_sub = C(surfacegreys_C+1:end,1);
%     for i = 1:size(C_info.brainstructurelabel,2)
%         thisstructsgreys = find(greybrainstructs_C ==i);
%         subcort_assings_by_structC{i,1} = C(thisstructsgreys);
%     end
% end
% %load('C:\Users\hermosir\Desktop\Test_data\MSC02a_to_ADHD315_template_MSC02b4_method_template_matching.mat')
% %eta_subject_index_4min = eta_subject_index;
% Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\Unique_5minute_ADHD_mode_paths.dscalar.nii',wb_command);
% Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\sub-33015a_task-rest_DCANBOLDProc_v4.0.0_Atlas_template_matched_Zscored.dscalar.nii',wb_command);
% %Dcii = ciftiopen('C:\Users\hermosir\Documents\test_ciftis\sub-100501_ses-20100430_avg_number_of_networks.dscalar.nii',wb_command);
%
% D = Dcii.cdata;
% if check_subcortical ==1
%     D_info = ft_read_cifti_mod('C:\Users\hermosir\Documents\test_ciftis\Unique_5minute_ADHD_mode_paths.dscalar.nii');
%     greybrainstructs_D = D_info.brainstructure(find(D_info.brainstructure >0));
%     surfacegreys_D = size(find(D_info.brainstructure <3 & D_info.brainstructure >0));
%     D_surf = D(1:surfacegreys_D,1);
%     D_sub = D(surfacegreys_D+1:end,1);
%     for i = 1:size(D_info.brainstructurelabel,2)
%         thisstructsgreys = find(greybrainstructs_D ==i);
%         subcort_assings_by_structD{i,1} = D(thisstructsgreys);
%     end
% end
% %test_missing_nets=1;
% if check_surface ==1
%     %D(find(D==10 | D==8 )) =5; % set all 10s to 1.
%     D=D_surf;C=C_surf;
% end



% numdat: number data. Each column is a catagory, each row represents
%   a separate set of data
% varargin{1}: cell of legend entries, one string for each column of numdat,
%   default is none, eg. {'First','Second','Third'}
% varargin{2}: cell of colors, one row of 3 RGB values for each category (column of numdat)
% varargin{3}: if 'pie', will make a standard pie chart
% varargin{4}: colormaska, will make a standard pie chart
% varargin{4}: colormaskb, will make a standard pie chart

% Examples:
%   donut([50,20,10;40,30,15],{'First','Second','Third'},{'r','g','b'});
%   donut([50,20,10],{'First','Second','Third'},[],'pie');
%   donut([50,20,10;40,30,15],[],{[ .945 .345 .329],[ .376 .741 .408],[ .365 .647 .855 ]});
% Default Values, if no variable arguments in
legtext = [];
colormap lines
clrmp = colormap;
ispie = 0;
% if length(varargin)>0
%     legtext = varargin{1};
%     if length(varargin)>1
%         if ~isempty(varargin{2})
%             clrmp = varargin{2};
%         else
%colormap lines

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
    200 200 200
    142 0 102]/255;

netRGBcs = netRGBs(isanet_in_current_setC,:);
netRGBds = netRGBs(isanet_in_current_setD,:);

if sort_mat ==1
    netRGBcs = netRGBcs(Eidx,:);
    netRGBds = netRGBds(Fidx,:);
else
end
clrmp = netRGBcs;
%clrmp = colormap;

%        end
%         if length(varargin)>2
%             if isempty(find(strcmp(varargin,'pie')))==0; ispie = 1; end
%         end
%                 if length(varargin)>2
%
%                 end
%     end
% end
rings = size(numdat,1); % nuber of rings in plot
%cats = size(numdat,2); % number of categories in each ring/set
%donout = nan(size(numdat{1},2));
for i = 1:rings
    for j = 1:size(numdat{i},2)
        percentage_vec(j) = (numdat{i}(j))/(sum(numdat{i}(:)))*100;
    end
    numdat_percentage{i,1} = percentage_vec;
end

for i = 1:rings
    tot = nansum(numdat{i,:}); % total things
    %donout(i,:)=numdat{i,:}./tot;
    fractang = (pi/2)+[0,cumsum((numdat{i,:}./tot).*(2*pi))];
    cats =  size(numdat{i,1},2);
    for j = 1:cats
        if ispie==1
            r0 = 0;
            r1 = 0.95;
        else
            r0 = i;
            r1 = i+0.95;
        end
        a0 = fractang(j);
        a1 = fractang(j+1);
        if iscell(clrmp)
            cl = clrmp{j};
        else
            if i ==1
                clrmp = netRGBcs;
                cl = clrmp(j,:);
                
            else
                clrmp = netRGBds;
                cl = clrmp(j,:);
            end
        end
        polsect(a0,a1,r0,r1,cl,skip_plotting);
        saved_polysect_values{i,j} = [a0,a1,r0,r1,cl];
    end
        set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[], 'YTick', []); 
    %set(gca,'Visible','off')
    axis off
    title(chart_title, 'Interpreter', 'none');
    
    %     if i==rings
    %         legend1 = legend(legtext);
    %         wi = legend1.Position(3);
    %         Xlm = xlim;
    %         widx = diff(Xlm);
    %         unitwi = widx.*wi;
    %         xlim([Xlm(1),Xlm(2)+unitwi])
    %     end
    
end

if add_label_percentages ==1
    for i = 1:rings
        for j = 1:cats
            try
            if numdat_percentage{i}(j)>label_threshold 
                %% new method to add percentages to the donuts
                %[a0,a1,r0,r1,c1,c2,c3] = saved_polysect_values{i,j};
                this_polsect = saved_polysect_values{i,j};
                a0 = this_polsect(1);
                a1 = this_polsect(2);
                r1 = this_polsect(4);
                
                
                
                labelRadius = 0.60 * r1; % reduce placement of label closer to the inside of the circle
                centerTheta = mean([a0 a1]);
                [halign,valign] = getAlignmentFromAngle(centerTheta);
                [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                
                % raw counts
                text(xtext,ytext,sprintf('%s',num2str(numdat_percentage{i}(j),3)),'fontsize',donut_label_font_size,'HorizontalAlignment',halign,'VerticalAlignment',valign);
            end
            catch
            end
        end
    end
end
end %function

function [halign, valign] = getAlignmentFromAngle(angle)
% Determine the text label alignment based on the angle around the circle.

% Convert the angle to degrees
angle = (180/pi)*angle;

% Round the angles to the nearest 45 degrees
angle = mod(round(angle/45)*45,360);

% Determine the horizontal alignment
if angle == 90 || angle == 270
    halign = 'center';
elseif angle > 90 && angle < 270
    halign = 'right';
else
    halign = 'left';
end

% Determine the vertical alignment
if angle == 0 || angle == 180
    valign = 'middle';
elseif angle > 180 && angle < 360
    valign = 'top';
else
    valign = 'bottom';
end

end


function pspatch = polsect(th0,th1,rh0,rh1,cl,skip_plotting)
% This function creates a patch from polar coordinates
a1 = linspace(th0,th0);
r1 = linspace(rh0,rh1);
a2 = linspace(th0,th1);
r2 = linspace(rh1,rh1);
a3 = linspace(th1,th1);
r3 = linspace(rh1,rh0);
a4 = linspace(th1,th0);
r4 = linspace(rh0,rh0);
[X,Y]=pol2cart([a1,a2,a3,a4],[r1,r2,r3,r4]);

if skip_plotting ==0
    
    p=patch(X,Y,cl); % Note: patch function takes text or matrix color def
    axis equal
    pspatch = p;
else
    pspatch =0;
end
end