function [all_rho, all_pval, all_nonmatched_num, all_nonmatched_percent] = make_percentage_correlation_plots(group1_conc_file,group2_conc_file,output_path,outputfilename_summary_name)

%add paths
this_code = which('make_percentage_correlation_plots');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);
disp('Attempting to add neccesaary paths and functions.')
warning('off') %supress addpath warnings to nonfolders.
for i=1:np
    addpath(genpath(settings.path{i}));
end
rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
rmpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti'); % remove non-working gifti path included with MSCcodebase
warning('on')
wb_command=settings.path_wb_c; %path to wb_command
group1_conc = importdata(group1_conc_file);
group2_conc = importdata(group2_conc_file);

for i=1:size(group1_conc)
    [~,B_temp,~] = fileparts(group1_conc{i});
    [~, outputfilename,~] = fileparts(B_temp);
    % changing network names
    network_match = regexp(outputfilename, '_(\w+)_network', 'tokens', 'once');
    if ~isempty(network_match)
        network_name = network_match{1};
        graph_network_cellname = strsplit(network_name,'_');
        graph_network_title = graph_network_cellname{end};
    else
        network_name = 'Unknown Network';
    end

    cii = ciftiopen(group1_conc{i},wb_command);
    dscalar_GRP1 = cii.cdata;
        
    cii = ciftiopen(group2_conc{i},wb_command);
    dscalar_GRP2 = cii.cdata;
    
    % check here to make sure that the dscalars that are loaded in are the same
    % length.  If one has a surface and one doesn't it will only compare
    % surfaces.
    if size(dscalar_GRP1,1) ~= size(dscalar_GRP2,1)
        if size(dscalar_GRP1,1) ==59412
            dscalar_GRP2 =dscalar_GRP2(1:59412,1);
        end
        if size(dscalar_GRP2,1) ==59412
            dscalar_GRP1 =dscalar_GRP1(1:59412,1);
        end
    end
    
    [rho_wzero,pval_wzero] = corr(double(dscalar_GRP1),double(dscalar_GRP2));
    
    %get indices of nonzero elements.
    log_nonZgrp1 = (dscalar_GRP1 ~=0);
    log_nonZgrp2 = (dscalar_GRP2 ~=0);
    
    both_nonzero = log_nonZgrp1 & log_nonZgrp2;
    
    x = (dscalar_GRP1(both_nonzero));
    y=(dscalar_GRP2(both_nonzero));
    
    c = ksdensity([x,y], [x,y]);
    
    % define figure properties
    opts.Colors     = get(groot,'defaultAxesColorOrder');
    opts.saveFolder = [output_path filesep];
    opts.width      = 8;
    opts.height     = 6;
    opts.fontType   = 'Times';
    opts.fontSize   = 20;
    
    
    figure();scatter(x, y, 10,c,'filled'); hold on
    plot([0:1],[0:1],'Color','k');
    
    set(gcf,'color','white')
    set(gca,'FontSize',20)
    caxis([ 1 20])
    cb = colorbar();
    colormap jet
    
    % fig.PaperPositionMode   = 'auto';
    title(graph_network_title);
    [rho,pval] = corr(x,y);
    performance = (y-x);
    performance_wzero = (dscalar_GRP1-dscalar_GRP2);
  
    mean_abs_error_wzero = mean(abs(performance_wzero)); 
    mean_abs_error = mean(abs(performance)); 
    
    log_Zerogrp1 = (dscalar_GRP1 ==0);
    log_Zerogrp2 = (dscalar_GRP2 ==0);
    disp(['The correlation between vectors is: ' num2str(rho) ', p=' num2str(pval) ])
    disp(['The correlation between vectors including zero elements is: ' num2str(rho_wzero) ', p=' num2str(pval_wzero) ])
    
    non_match_zero = xor(log_Zerogrp1,log_Zerogrp2);
    disp(['Number of greyordinates with mismatched zeros is: ' num2str(sum(non_match_zero)) ' (or ' num2str((sum(non_match_zero))/size(both_nonzero,1)) '% of the greyordinates).'])
    
    disp(['Saving image: ' opts.saveFolder outputfilename])
    print([opts.saveFolder outputfilename], '-dpng', '-r600')
    
    
    all_rho(i) = rho;
    all_pval(i) = pval;
    all_rho_wzero(i) = rho_wzero;
    all_pval_wzero(i) = pval_wzero;
    all_performance{i} = performance;
    all_performance_wzero{i} = performance_wzero;
    all_mean_abs_error(i) = mean_abs_error;
    all_mean_abs_error_wzero(i) = mean_abs_error_wzero;
   
    all_nonmatched_num(i) = sum(non_match_zero);
    all_nonmatched_percent(i) = (sum(non_match_zero))/size(both_nonzero,1);
    cii.cdata = performance_wzero;
    ciftisave(cii,[output_path filesep outputfilename '_performance.dscalar.nii'],wb_command);
end
disp('Saving .mat file with summary stats saving...')
save([output_path filesep outputfilename '_correlation_.mat'],'all_rho','all_pval','all_nonmatched_num','all_nonmatched_percent','all_rho_wzero','all_pval_wzero','all_performance','all_performance_wzero','all_mean_abs_error','all_mean_abs_error_wzero')
disp('Done running code for all files in conc.')

end
