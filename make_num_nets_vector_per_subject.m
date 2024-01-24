function [whole_brain_number_of_nets_all,integration_zone_number_of_nets_all] = make_num_nets_vector_per_subject(dtseriesconc, outputdir,output_name)
%This function generates a vector of the number of networks for each
%subject based on a dtseries with overlapping network labels.


dtserieswithassignments = importdata(dtseriesconc);

%preallocate space for saved variables.
integration_zone_number_of_nets_all = zeros(length(dtserieswithassignments),1);
whole_brain_number_of_nets_all = zeros(length(dtserieswithassignments),1);
all_uniqueIDs = cell(length(dtserieswithassignments),1);

tic
for i = 1:length(dtserieswithassignments)
    disp(i)
    if rem(i,100)==0 % provide time updates
        disp([' Processing time update: Subject: ' num2str(i)]);toc;
    end
    dtseriesname = dtserieswithassignments{i};
    %get output_names
    %HARDCODE WARNING
    %uniqueID = char(dtseriesname(88:110));
    uniqueID = char(dtseriesname(108:126)); % ABCD MSI
    %uniqueID = char(dtseriesname(101:123));
    %uniqueID = char(dtseriesname(96:113)); %HCP-D dtseries
    %uniqueID = char(dtseriesname(99:108)); %BCP8-11
    %uniqueID = char(dtseriesname(100:109)); %BCP8-11
    
    disp(uniqueID)
    all_uniqueIDs{i} = uniqueID;
    outputname = [outputdir filesep uniqueID];
    if exist([outputname '_avg_number_of_networks.dscalar.nii'],'file') ~=0
        disp([outputname '_avg_number_of_networks.dscalar.nii found. loading .mat data.'])
        load([outputname '.mat'],'whole_brain_number_of_nets','integration_zone_number_of_nets');
    else
        [~,whole_brain_number_of_nets,integration_zone_number_of_nets] = visualizedscalars(dtserieswithassignments{i},outputname,'number_of_networks',0,0,0);
    end
    whole_brain_number_of_nets_all(i) = whole_brain_number_of_nets;
    integration_zone_number_of_nets_all(i) = integration_zone_number_of_nets;
    
    clear whole_brain_number_of_nets integration_zone_number_of_nets % clear just in case
end
save([output_name '.mat'],'whole_brain_number_of_nets_all','integration_zone_number_of_nets_all','dtserieswithassignments','dtseriesconc','all_uniqueIDs');

disp('Done making number of networks vectors per subject');

