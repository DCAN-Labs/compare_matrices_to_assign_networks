function [outname] = clean_dscalars_by_size(dscalarwithassignments,manualset,groupnetworksfile,dostripes,mincol,minsize,orig_parcelsfile,make_consensus,assign_unassigned,remove46)
%consensus_maker_knowncolors(regularized_ciftifile,[manualset],[groupnetworksfile],[dostripes],[mincol],[minsize],[orig_parcelsfile])

%This function cleans up networks.
%Hardcode
%make_consensus = 1;

conc = strsplit(dscalarwithassignments, '.');
conc = char(conc(end));
if strcmp('conc',conc) == 1
    dscalarwithassignments = importdata(dscalarwithassignments);
else
    dscalarwithassignments = {dscalarwithassignments};
end

%check to make sure that surface files exist
tic
for i = 1:length(dscalarwithassignments)
    if rem(i,100)==0
        disp([' Validating file existence ' num2str(i)]);toc;
    end
    if exist(dscalarwithassignments{i}) == 0
        NOTE = ['Subject dscalar ' num2str(i) ' does not exist']
        disp(dscalarwithassignments{i});
        return
    else
    end
end
disp('All series files exist continuing ...')

network_assignment_filetype = strsplit(dscalarwithassignments{i}, '.');
cifti_type = char(network_assignment_filetype(end-1));
if strcmp('dtseries',cifti_type) == 1
    overlap =1;
else
    overlap =0;
end

%% Adding paths for this function
this_code = which('clean_dscalars_by_size');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
%support_folder=[pwd '/support_files'];
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

if isdeployed
    disp('Matlab is deployed. Not adding paths, as they should have been added during compiling.')
else
    warning('off') %supress addpath warnings to nonfolders.
    for i=1:np
        addpath(genpath(settings.path{i}));
    end
    rmpath('/mnt/max/shared/code/external/utilities/MSCcodebase/Utilities/read_write_cifti') % remove non-working gifti path included with MSCcodebase
    rmpath('/home/exacloud/fnl_lab/code/external/utilities/MSCcodebase/Utilities/read_write_cifti')
    addpath(genpath('/home/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/'));
    wb_command=settings.path_wb_c; %path to wb_command
    warning('on')
end

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('dostripes','var') || isempty(dostripes)
    dostripes = 0;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end

if ~exist('groupnetworksfile','var') || isempty(groupnetworksfile)
    groupnetworksfile = settings.path{6}; %Networks_template_cleaned.dscalar.nii
    %groupnetworksfile = '/home/data/atlases/Networks_template.dscalar.nii';
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:100];


%% Start code
tic
for i=1:length(dscalarwithassignments)
    if rem(i,100)==0
        disp([' Cleaning subject ' num2str(i)]);toc;
    end
    regularized_ciftifile = dscalarwithassignments {i};
    
    cifti_data = ft_read_cifti_mod(regularized_ciftifile); assigns = cifti_data.data;

    assigns(assigns<0) = 0; %Assignments that are eqaul to -1 (unassigned), set the to 0.
    assigns(isnan(assigns)) = 0; % Set nans to 0.
    
    if make_consensus == 1
        groupfile = ft_read_cifti_mod(groupnetworksfile);
        groupdata = groupfile.data;
        ncortverts = nnz(groupfile.brainstructure==1) + nnz(groupfile.brainstructure==2); %grab the cortical vetices of the template for jaccard.
        groupdata = groupdata(1:ncortverts,1);
        
        potential_colors = [1 2 10 9 3 5 11 16 15 7 8 12 14 13];%[1 2 10 9 3 5 6 11 16 15 7 8 17 12 4 14 13];
        newcolors = setdiff(all_color_values,potential_colors); % these are the potential network assignsment for "non-canonical " networks.
        
        unassigned_networks = cell(1,size(assigns,2));
        
        all_recolored = zeros(size(assigns)); % preallocate columns
        for c = 1:size(all_recolored,2) 
            col_consensusmap = assigns(:,c);
            
            
            unassigned = find(col_consensusmap<1);
            for unassignedindex = unassigned'
                thisassignments = assigns(unassignedindex,mincol:end);
                thisassignments(thisassignments<1) = [];
                if ~isempty(thisassignments)
                    col_consensusmap(unassignedindex) = thisassignments(1);
                end
            end
            networks = unique(col_consensusmap); networks(networks<=0) = [];
            new_networks = networks;
            assigning_networks = networks;
            
            col_out = zeros(size(col_consensusmap));
            
            if exist('manualset')
                if isempty(manualset)
                    clear manualset
                else
                    for i = 1:size(manualset,1)
                        col_out(col_consensusmap==manualset(i,1)) = manualset(i,2);
                        new_networks(new_networks==manualset(i,1)) = [];
                        if all(manualset(i,2)~=potential_colors); new_networks = [new_networks; manualset(i,2)]; end
                        assigning_networks(assigning_networks==manualset(i,1)) = [];
                    end
                end
            end
            
            col_consensusmap_nosubcort = col_consensusmap(1:size(groupdata,1));
            
            for i = 1:length(potential_colors)
                
                if exist('manualset') && any(manualset(:,2)==potential_colors(i))
                    
                    
                    
                else
                    
                    
                    if ~isempty(assigning_networks)
                        groupnetwork_comp = groupdata==potential_colors(i);
                        D = zeros(length(assigning_networks),1);
                        P = zeros(length(assigning_networks),1);
                        for j = 1:length(assigning_networks)
                            
                            network_comp = col_consensusmap_nosubcort==assigning_networks(j);
                            P(j) = nnz(groupnetwork_comp & network_comp);
                            D(j) = P(j)/nnz(groupnetwork_comp | network_comp);
                            
                            if c>2 && any(any(all_recolored(1:ncortverts,1:(c-1))==potential_colors(i)))
                                prevnetwork_comp = any(all_recolored(1:ncortverts,1:(c-1))==potential_colors(i),2);
                                D_withprev = nnz(prevnetwork_comp & network_comp) ./ nnz(network_comp | prevnetwork_comp);
                                if D_withprev < .1
                                    D(j) = 0;
                                end
                            end
                        end
                        [maxval, maxind(i)] = max(D);
                        
                        if maxval > .1
                            col_out(col_consensusmap==assigning_networks(maxind(i))) = potential_colors(i);
                            new_networks(new_networks==assigning_networks(maxind(i))) = [];
                            assigning_networks(assigning_networks==assigning_networks(maxind(i))) = [];
                        end
                    end
                end
                
            end
            clear maxind D P
            for j = 1:length(new_networks)
                col_out(col_consensusmap==new_networks(j)) = newcolors(j);
            end
            all_recolored(:,c) = col_out;
            unassigned_networks{c} = assigning_networks;
        end
        
        
        
        
        
        all_recolored(assigns<=0) = 0;
        
        if remove46 == 1 % set values that are either 4 or 6 to zero (for missing networks).
            
            all_recolored(all_recolored==4) = 0;
            all_recolored(all_recolored==6) = 0;
            
        end
        
        
        cifti_data.data = all_recolored;
        if ~exist('cifti_data.mapname')
            for col = 1:size(all_recolored,2)
                cifti_data.mapname{col} = ['Column number ' num2str(col)];
            end
            cifti_data.dimord = 'scalar_pos';
        end
        
        
        
        dotsloc = strfind(regularized_ciftifile,'.');
        basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
        outname = [basename '_allcolumns_recolored'];
        ft_write_cifti_mod(outname,cifti_data);
        set_cifti_powercolors([outname '.dscalar.nii'])
        
        
        
        
        out = all_recolored(:,mincol);
        
        uniquevals = unique(out); %uniquevals(uniquevals<1) = [];
        colors_tofix = setdiff(uniquevals,potential_colors);
        verts_tofix = [];
        for colornum = 1:length(colors_tofix)
            verts_thiscolor = find(out==colors_tofix(colornum));
            verts_tofix = [verts_tofix ; verts_thiscolor(:)];
        end
        for vertnum = verts_tofix'
            for col = (mincol+1):size(all_recolored,2)
                if any(all_recolored(vertnum,col)==potential_colors)
                    out(vertnum) = all_recolored(vertnum,col);
                    break
                end
            end
        end
    else
        out = assigns;
    end
    temp_out = out;
    temp_all =  temp_out;
    if logical(minsize)
        
        
        % Clean up tiny twinspieces
        if size(temp_out,1) == 91282
            if exist([support_folder '/node_neighbors_91282.mat'],'file') == 2
                load([support_folder '/node_neighbors_91282.mat']);
            else
                %neighbors = cifti_neighbors(regularized_ciftifile); Luci
                %change 6/10/2021, replaced with following line:
                neighbors = cifti_neighbors_dcan(regularized_ciftifile,[],[],path_to_neighbors_file);
            end
        else
             if exist([support_folder '/node_neighbors_59412.mat'],'file') == 2
              load([support_folder '/node_neighbors_59412.mat']);
            else
                %error([char([support_folder '/node_neighbors_59412.mat']) ' file is missing.'])
                %neighbors = cifti_neighbors(regularized_ciftifile); Luci
                %change 6/10/2021: removed this line and replaced with
                %following 2 lines in order to work for surface only data:
                %path_to_neighbors_file = [support_folder filesep 'node_neighbors.txt'];
                path_to_neighbors_file = [support_folder filesep 'node_neighbors_59412.mat']; 
                neighbors = cifti_neighbors_dcan(regularized_ciftifile,[],[],path_to_neighbors_file);
            end
        end
        
        for j = 1:size(out,2)
            if j==4 || j ==6
                continue
            end
            
            allcolors= unique(out(:,j));
            if assign_unassigned == 1
                %allcolors = [allcolors(2:end)' allcolors(1)']'; % move unsigned networks to the end, so they're likely to be already assinged.
            else %leave unsigned vertices unassinged.
                allcolors(allcolors<=0) = [];
            end
            temp_out = temp_all(:,j);
            
            
            for color = allcolors(:)'
                clusteredmetric = zeros(size(temp_out));
                thiscolorverts = find(temp_out==color);
                for vertex = thiscolorverts'
                    
                    %find the neighbors of this vertex
                    vertexneighbors = neighbors(vertex,:);
                    vertexneighbors(isnan(vertexneighbors)) = [];
                    
                    %find which of those neighbors also pass the thresholds
                    vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
                    
                    %find if those neighbors have already been assigned different cluster values
                    uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
                    
                    uniqueneighborvals(uniqueneighborvals==0) = [];
                    
                    
                    %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
                    if isempty(uniqueneighborvals)
                        clusteredmetric(vertexneighbors_thiscolor) = vertex;
                        %if there is only one previous cluster identifier present, make all the neighbors that value
                    elseif length(uniqueneighborvals)==1
                        clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
                        %if there are multiple cluster identifier values in the neighborhood, merge them into one
                    else
                        for valuenum = 2:length(uniqueneighborvals)
                            clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                            if color == 0
                                if assign_unassigned ==1 % added by Robert H.
                                    clusteredmetric(vertex) = uniqueneighborvals(1); % don't forget to assign the vertex the same assingment.
                                else
                                end
                            else
                            end
                            
                        end
                    end
                end
                uniqueclustervals = unique(clusteredmetric);
                uniqueclustervals(uniqueclustervals==0) = [];
                
                
                
                for clusternum = uniqueclustervals' %added by Robert
                    %                disp(['cluster ID ' num2str(clusternum)])
                    verts_in_cluster = nnz(clusteredmetric==clusternum);
                    %                     disp([num2str(verts_in_cluster) ' vertices in cluster'])
                    %                end
                    
                    %if nnz(clusteredmetric==clusternum) < minsize
                    if verts_in_cluster < minsize
                        neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                        neighborverts(isnan(neighborverts)) = [];
                        borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                        if size(temp_out,1) ~= 91282
                        large_val=find(borderverts>59412);
                        borderverts(large_val)=NaN;
                        borderverts=borderverts(~isnan(borderverts));
                        end
                        %borderverts(temp_out(borderverts)<1) = [];
                        %added by robert
                        if assign_unassigned ==1
                            %borderverts(temp_out(borderverts)) = [];
                            bordererassigns = temp_out(borderverts);
                            int_clusterassings = bordererassigns > 0;
                            mode_neighborval = mode(bordererassigns(int_clusterassings));
                            if mode_neighborval ==0
                                disp(clusternum)
                            end
                        else
                            if overlap ==0 % allow networks to get an assingment of 0, but only for overlapping networks.
                                borderverts(temp_out(borderverts)<1) = [];
                                %mode_neighborval = mode(temp_out(isassinged(borderverts)));
                            else
                                borderverts(temp_out(borderverts)<0) = [];
                            end
                            mode_neighborval = mode(temp_out(borderverts));
                            %                             if isnan(mode_neighborval) ==1
                            %                                 mode_neighborval =0;
                            %                             end
                        end
                        %Grab the next value
                        %mode_neighborval = mode(temp_out(borderverts));
                        temp_out(clusteredmetric==clusternum) = mode_neighborval;
                    else % don't forget to assign large clusters to networks (network ==0) to values is still unassigned.
                        if assign_unassigned ==1
                            if color == 0
                                neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                                neighborverts(isnan(neighborverts)) = [];
                                borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                                bordererassigns = temp_out(borderverts);
                                int_clusterassings = bordererassigns > 0;
                                mode_neighborval = mode(bordererassigns(int_clusterassings));
                                temp_out(clusteredmetric==clusternum) = mode_neighborval;
                                
                            else
                            end
                        else
                        end
                    end
                end
                
                
            end
            all_temp_out(:,j) = temp_out;
            
        end
        out = temp_out;
        
    end
    
    if overlap == 1
        cifti_data.data = all_temp_out;
    else
        cifti_data.data = out;
    end
    
    
    if overlap == 0
        if ~exist('cifti_data.mapname')
            cifti_data.mapname = {'Column number ' num2str(mincol)};
            cifti_data.dimord = 'scalar_pos';
        else
            cifti_data.mapname = cifti_data.mapname(mincol);
        end
    else
        if ~exist('cifti_data.mapname')
            cifti_data.mapname = {'Column number ' num2str(mincol)};
            cifti_data.dimord = 'pos_time';
        else
            cifti_data.mapname = cifti_data.mapname(mincol);
        end
    end
    dotsloc = strfind(regularized_ciftifile,'.');
    basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
    outname = [basename '_recolored'];
    ft_write_cifti_mod(outname,cifti_data);
    if overlap ==0
        
        set_cifti_powercolors([outname '.dscalar.nii'])
        
    else
        set_cifti_powercolors([outname '.dtseries.nii'])
        
    end
    if exist('orig_parcelsfile') && ~isempty(orig_parcelsfile)
        parcels = ft_read_cifti_mod(orig_parcelsfile);
        parcels = parcels.data;
        parcels((length(out)+1):end) = [];
        IDs = unique(parcels); IDs(IDs<1) = [];
        outtext = zeros(length(IDs),1);
        outtext_bycol = zeros(length(IDs),size(all_recolored,2));
        for IDnum = 1:length(IDs)
            outtext(IDnum) = mode(out(parcels==IDs(IDnum)));
            outtext_bycol(IDnum,:) = mean(all_recolored(parcels==IDs(IDnum),:),1);
        end
        dlmwrite([outname '.txt'],outtext,'delimiter',' ')
        dlmwrite([basename '_allcolumns_recolored.txt'],outtext_bycol,'delimiter',' ')
    end
end
disp('Done running clean dscalars code.')

%Stripes

% if dostripes
%
%     all_recolored = all_recolored(1:ncortverts,:);
%     allcolors = unique(all_recolored); allcolors(allcolors==0) = [];
%     unknown_colors = setdiff(allcolors,potential_colors);
%     for color = unknown_colors(:)'
%         all_recolored(all_recolored==color) = 0;
%     end
%
%     to_be_striped = out;
%     change = logical(diff(all_recolored,1,2));
%     change = change .* (all_recolored(:,1:end-1)>0) .* (all_recolored(:,2:end)>0);
%     for col = 1:size(change,2)
%         colvals = unique(all_recolored(:,col+1)); colvals(colvals==0) = [];
%         for val_totest = colvals(:)'
%             if nnz((all_recolored(:,col)==val_totest) & (all_recolored(:,col+1)==val_totest)) / nnz((all_recolored(:,col)==val_totest) | (all_recolored(:,col+1)==val_totest)) > .5
%                 thiscolor_changed_inds = logical(change(:,col) .* (all_recolored(:,col)==val_totest));
%                 to_be_striped(thiscolor_changed_inds,2+((col-1)*2)) = all_recolored(thiscolor_changed_inds,col+1);
%                 to_be_striped(thiscolor_changed_inds,3+((col-1)*2)) = val_totest;
%             end
%         end
%     end
%
%
%     to_be_striped_final = zeros(size(to_be_striped));
%     for vert = 1:size(to_be_striped_final,1)
%         uniquevals = unique([out(vert) to_be_striped(vert,:)]); uniquevals(uniquevals==0) = [];
%         to_be_striped_final(vert,1:length(uniquevals)) = uniquevals;
%     end
%     make_striped_cifti(to_be_striped_final,0,[outname '_striped_164'],1/30)
%     set_cifti_powercolors([outname '_striped_164.dtseries.nii'])
%
% end

