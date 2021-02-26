function [outputname] = concatenateFDnumbers(FD_conc_file,outputname)
%This function works by loading in a list of FD values (usually listed as
%motion_numbers.mat) and concatnating them.

% if exist('FD_conc_file', 'var') == 1 && ~isempty(FD_conc_file) == 1
%     if strcmp('conc',conc) == 1
         FD_conc_file = importdata(FD_conc_file);
%     else
%         FD_conc_file = {FD_conc_file};
%     end
% else
% end

for i = 1:length(FD_conc_file)
    load(FD_conc_file{i})
    
    current_diff_x = motion_numbers.diff_x;
    current_diff_y = motion_numbers.diff_y;
    current_diff_z = motion_numbers.diff_z;
    current_diff_arc_len_x = motion_numbers.diff_arc_len_x;
    current_diff_arc_len_y = motion_numbers.diff_arc_len_y;
    current_diff_arc_len_z = motion_numbers.diff_arc_len_z;
    current_FD = motion_numbers.FD;
    current_DVAR_pre_reg = motion_numbers.DVAR_pre_reg;
    current_DVAR_post_reg = motion_numbers.DVAR_post_reg;
    current_DVAR_post_all = motion_numbers.DVAR_post_all;
    %     mask = importdata(all_motion_conc{i});
    %     FD_afterscrub = FD(logical(mask));
    %     avgFD_prescrub(i,:) =mean(FD);
    %     avgFD_postscrub(i,:) = mean(FD_afterscrub);
    if i ==1
        
        diff_x_all =  current_diff_x;
        diff_y_all = current_diff_y;
        diff_z_all = current_diff_z;
        diff_arc_len_x_all = current_diff_arc_len_x;
        diff_arc_len_y_all = current_diff_arc_len_y;
        diff_arc_len_z_all = current_diff_arc_len_z;
        FD_all  = current_FD;
        DVAR_pre_reg_all = current_DVAR_pre_reg;
        DVAR_post_reg_all = current_DVAR_post_reg;
        DVAR_post_all_all = current_DVAR_post_all;
        
    else
        %concatenate values
        diff_x_all = [diff_x_all; current_diff_x];
        diff_y_all = [diff_y_all; current_diff_y];
        diff_z_all = [diff_z_all; current_diff_z];
        diff_arc_len_x_all = [diff_arc_len_x_all; current_diff_arc_len_x];
        diff_arc_len_y_all = [diff_arc_len_y_all; current_diff_arc_len_y];
        diff_arc_len_z_all = [diff_arc_len_z_all; current_diff_arc_len_z];
        FD_all = [FD_all; current_FD];
        DVAR_pre_reg_all = [DVAR_pre_reg_all; current_DVAR_pre_reg];
        DVAR_post_reg_all = [DVAR_post_reg_all; current_DVAR_post_reg];
        DVAR_post_all_all = [DVAR_post_all_all; current_DVAR_post_all];
        
    end   
end

%write all the files
    motion_numbers.diff_x = diff_x_all;
    motion_numbers.diff_y = diff_y_all;
    motion_numbers.diff_z = diff_z_all;
    motion_numbers.diff_arc_len_x = diff_arc_len_x_all;
    motion_numbers.diff_arc_len_y = diff_arc_len_y_all;
    motion_numbers.diff_arc_len_z = diff_arc_len_z_all;
    motion_numbers.FD = FD_all;
    motion_numbers.DVAR_pre_reg = DVAR_pre_reg_all;
    motion_numbers.DVAR_post_reg = DVAR_post_reg_all;
    motion_numbers.DVAR_post_all = DVAR_post_all_all;

disp(['saving motion numbers to: ' outputname])
save(outputname,'motion_numbers')


