function h = alluvialflow_change_colors(data, left_labels, right_labels, chart_title,sort_mat, Eidx, Fidx,isanet_in_current_setC, isanet_in_current_setD,output_name)


%
% Plot an alluvial flow diagram.
% left_labels:  Names of categories to flow from.
% right_labels: Names of categories to flow to.
% data:         Matrix with size numel(left_labels) rows by
%               numel(right_labels) columns.
%
% Ideas for future work:
% 1. Get data from a MATLAB table, use table variable names and named rows
%    as labels.
% 2. Interface similar to the plot function, with name-value pairs, optional
%    parameters etc.
h = gcf;
%clf
% set(h, 'WindowStyle', 'Docked'); % DBG this helps reuse desktop space
set(gcf,'color','w');
% Find axis dimensions and set them
data_sum = sum(data(:));
total_gap = 0.10 * data_sum;

left_gap_size = total_gap / (size(data, 1)-1);

right_gap_size = total_gap / (size(data, 2)-1);
if size(data,1) ==1
    left_gap_size = 0;
end
if size(data,2) ==1
    right_gap_size = 0;
end
y_height = data_sum + total_gap;
x_left = 0;
x_right = 1;
axis([x_left, x_right, 0, y_height]) % Set limits
axis ij % origin is top left
axis off

% grid minor % DBG

hold on
patch([0 0 1 1], [0 y_height y_height 0], 'w');

% Plot left categories - one per row
left_category_sizes = sum(data, 2)';

% These are the top points for each left category,
% with gaps added.
left_category_points = [0 cumsum(left_category_sizes)] + ...
    (0:numel(left_category_sizes)) .* left_gap_size;
left_category_points(end) = [];

% plot left category bars
plot(zeros(2, numel(left_category_points)), [left_category_points; (left_category_points + left_category_sizes)], 'k', 'LineWidth',4);

% DBG plot left ticks
%left_category_tick_starts = zeros(size(left_category_points)) - 0.01;
%left_category_tick_ends = left_category_tick_starts + 0.02;
%plot([left_category_tick_starts; left_category_tick_ends], ...
%     [left_category_points; left_category_points], 'b-');

% Plot right categories - one per column
right_category_sizes = sum(data, 1);

% These are the top points for each right category,
% with gaps added.
right_category_points = [0 cumsum(right_category_sizes)] + ...
    (0:numel(right_category_sizes)) .* right_gap_size;
right_category_points(end) = [];

% plot right category bars
plot(ones(2, numel(right_category_points)), [right_category_points; (right_category_points + right_category_sizes)], 'k', 'LineWidth',4);

% DBG plot right ticks
%right_category_tick_ends = ones(size(right_category_points)) + 0.01;
%right_category_tick_starts = right_category_tick_ends - 0.02;
%plot([right_category_tick_starts; right_category_tick_ends], ...
%     [right_category_points; right_category_points], 'b-');

%
% Draw the patches, an entire left category at a time
%

% Color selection
patch_colors = [ .5 .5 .5;
    1  0  0;
    0  1  0;
    0  0  1;
    .5 .5  0;
    0 .5 .5;
    .5  0 .5];
num_colors = size(patch_colors, 1);

right_columns_so_far = right_category_points(1:end); % Start at the beginning of each right category and stack as we go.
patches_per_left_category = size(data, 2);

for k_left = 1:size(data, 1) % for each row
    color = patch_colors(mod(k_left,num_colors)+1, :);
    
    %
    % Calculate the coordinates for all the patches split by the
    % current left category
    %
    
    % Split the left category
    left_patch_points = [0 cumsum(data(k_left, :))] + left_category_points(k_left);
    patch_top_lefts = left_patch_points(1:end-1);
    patch_bottom_lefts = left_patch_points(2:end);
    
    % Compute and stack up slice of each right category
    patch_top_rights = right_columns_so_far;
    patch_bottom_rights = patch_top_rights + data(k_left, :);
    right_columns_so_far = patch_bottom_rights;
    
    %
    % Plot the patches
    %
    
    % X coordinates of patch corners
    [bottom_curves_x, bottom_curves_y] = get_curves(0.1, patch_bottom_lefts, 0.9, patch_bottom_rights);
    [top_curves_x,    top_curves_y]    = get_curves(0.9, patch_top_rights,   0.1, patch_top_lefts);
    X = [ ...
        repmat([0; 0], 1, patches_per_left_category); % Top left, bottom left
        bottom_curves_x;
        repmat([1; 1], 1, patches_per_left_category); % Bottom right, top right
        top_curves_x
        ];
    %X1col = X(:,1);
    %X2col = X(:,2);
    % Y coordinates of patch corners
    Y = [ ...
        patch_top_lefts;
        patch_bottom_lefts;
        bottom_curves_y;
        patch_bottom_rights;
        patch_top_rights;
        top_curves_y
        ];
    %Y1col = Y(:,1);
    %Y2col = Y(:,2);
    %test_colors =mymap(1:size(X),1:3);
    %     load('PowerColorMap.mat')
    
    %     for RGB_columns = 1:3
    %         for c = 1:(size(mymap,1)):8
    %             net_colors(:,RGB_columns) = linspace(mymap(c,RGB_columns),mymap((c+7),RGB_columns),size(X,1));
    %         end
    %     end
%    if overlap ==0
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
        142 0 103]/255;
%     else
%         netRGBs = [
%         255 0 0;
%         0 0 153
%         255 255 0
%         0 255 0
%         13 133 160
%         50 50 50
%         102 0 204
%         102 255 255
%         255 128 0
%         178 102 153
%         0 102 153
%         102 255 102
%         60 60 251
%         200 200 200
%         0 0 0 ]/255;  
%     end
    netRGBcs = netRGBs(isanet_in_current_setC,:);
    netRGBds = netRGBs(isanet_in_current_setD,:);
    
    if sort_mat ==1
        netRGBcs = netRGBcs(Eidx,:);
        netRGBds = netRGBds(Fidx,:);
    else
    end
    %          net_colors = [
    %             repmat(2,1,size(patch_top_lefts,2));
    %             repmat(2,1,size(patch_bottom_lefts,2));
    %             repmat(0,15,size(bottom_curves_y,2));
    %             repmat(3,1,size(patch_bottom_rights,2));
    %             repmat(3,1,size(patch_top_rights,2));
    %             repmat(0,15,size(top_curves_y,2));
    %              ];
    for L=1:size(X,2)
        
        for RGB_columns = 1:3
            %bottom_curve_RGBS(:,RGB_columns)  = linspace(netRGBs(1,RGB_columns),netRGBs(2,RGB_columns),size(bottom_curves_x,1));
            %if sort_mat ==0
            %bottom_curve_RGBS(:,RGB_columns)  = linspace(netRGBs(k_left,RGB_columns),netRGBs(L,RGB_columns),size(bottom_curves_x,1));
            %else
            bottom_curve_RGBS(:,RGB_columns)  = linspace(netRGBcs(k_left,RGB_columns),netRGBds(L,RGB_columns),size(bottom_curves_x,1));
            %end
        end
            
            for RGB_columns = 1:3
                %top_curve_RGBS(:,RGB_columns)  = linspace(netRGBs(2,RGB_columns),netRGBs(1,RGB_columns),size(bottom_curves_x,1));
                %if sort_mat ==0
                %top_curve_RGBS(:,RGB_columns)  = linspace(netRGBs(L,RGB_columns),netRGBs(k_left,RGB_columns),size(bottom_curves_x,1));
                %else
                top_curve_RGBS(:,RGB_columns)  = linspace(netRGBds(L,RGB_columns),netRGBcs(k_left,RGB_columns),size(bottom_curves_x,1));
                %end
            end
        
       if sort_mat ==1
        
        net_per_flow= [netRGBcs(k_left,:);
            %net_per_flow= [netRGBs(k_left,:);
            
            %netRGBs(k_left,:);
            netRGBcs(k_left,:);
            
            bottom_curve_RGBS
            netRGBds(L,:);
            %netRGBs(L,:);
            netRGBds(L,:);
            %netRGBs(L,:);
            top_curve_RGBS;
            %bottom_curve_RGBS
            ];
        
       else
          net_per_flow= [netRGBcs(k_left,:);
            
            netRGBcs(k_left,:);
            %netRGBas(k_left,:);
            
            bottom_curve_RGBS
            %netRGBbs(L,:);
            netRGBds(L,:);
            %netRGBbs(L,:);
            netRGBds(L,:);
            top_curve_RGBS;
            %bottom_curve_RGBS
            ];
        
       end
        
        one_flow_colors=zeros(34,1,3);
        one_flow_colors(:,1,:)= net_per_flow;
        
        %g= patch( X2col, Y2col,one_flow_colors,'FaceAlpha', .7, 'EdgeColor', 'none');
        
        
        patch( X(:,L), Y(:,L),one_flow_colors,'FaceAlpha', .7, 'EdgeColor', 'none');
    end
    %g= patch('XData', X, 'YData', Y, 'FaceColor',patch_colors(k_left,:), 'FaceAlpha', .4, 'EdgeColor', 'none');
    %g= patch('XData', X, 'YData', Y, 'FaceColor',test_colors(1,1:3), 'FaceAlpha', .4, 'EdgeColor', 'none');
    %colormap(mymap)
    
    
    
end % for each row

% Place left labels
text(zeros(1, size(data, 1)) - 0.01, ...
    left_category_points + left_category_sizes./2, ...
    left_labels, 'FontSize', 7, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);

% Place right labels
text(ones(1, size(data, 2)) + 0.01, ...
    right_category_points + right_category_sizes./2, ...
    right_labels, 'FontSize', 7, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Rotation', 90);

title(chart_title, 'Interpreter', 'none');
h.Position = [100 100 200 700]; %make a 100 x 600 figure

 print('-dpng','-r300',h,[output_name '_alluvial_fig'])

end % alluvialflow

function [x, y] = get_curves(x1, y1, x2, y2)
% x1, x2: scalar x coordinates of line start, end
% y1, y2: vectors of y coordinates of line start/ends
Npoints = 15;
t = linspace(0, pi, Npoints);
c = (1-cos(t))./2; % Normalized curve

Ncurves = numel(y1);
% Starting R2016b, the following line could be written simply as:
%   y = y1 + (y2 - y1) .* c';
y = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints,1) .* repmat(c', 1, Ncurves);
x = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end  % get_curve