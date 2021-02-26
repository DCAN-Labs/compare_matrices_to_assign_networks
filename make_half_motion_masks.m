
function make_half_motion_masks(outputfilename,frames)

% This function makes "additional masks" ( i.e. a vector 1s and 0s). to use
% with the cifti_conn_matrix function.

maska = ones(frames,1);
maskb = ones(frames,1);

maska(round((frames/2))+1:end) = 0;
maskb(round(1:frames/2)) = 0; 

fileID = fopen([outputfilename '_additionalmaska.txt'],'w');
fprintf(fileID,'%d\n',maska);
fclose(fileID);

fileID = fopen([outputfilename '_additionalmaskb.txt'],'w');
fprintf(fileID,'%d\n',maskb);
fclose(fileID);

end