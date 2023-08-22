
function make_half_motion_masks(outputfilename,num_frames)

% This function makes "additional masks" ( i.e. a vector 1s and 0s). to use
% with the cifti_conn_matrix function.

maska = ones(num_frames,1);
maskb = ones(num_frames,1);

maska(round((num_frames/2))+1:end) = 0;
maskb(round(1:num_frames/2)) = 0; 

fileID = fopen([outputfilename '_additionalmaska.txt'],'w');
fprintf(fileID,'%d\n',maska);
fclose(fileID);

fileID = fopen([outputfilename '_additionalmaskb.txt'],'w');
fprintf(fileID,'%d\n',maskb);
fclose(fileID);

end