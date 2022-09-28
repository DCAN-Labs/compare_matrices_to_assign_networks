function [whole_brain_number_of_nets_all,integration_zone_number_of_nets_all] = make_num_nets_vector_per_subject(dtseriesconc, outputdir,output_name)
%This function generates a vector of the number of networks for each
%subject based on a dtseries with overlapping network labels.

dtserieswithassignments = importdata(dtseriesconc);
integration_zone_number_of_nets_all = zeros(size(length(dtserieswithassignments),1));
whole_brain_number_of_nets_all = zeros(size(length(dtserieswithassignments),1));

for i = 1:length(dtserieswithassignments)
    disp(i)
    dtseriesname = dtserieswithassignments{i};
    %get output_names
    %HARDCODE WARNING
    %uniqueID = char(dtseriesname(88:110));
    %uniqueID = char(dtseriesname(101:123));
    %uniqueID = char(dtseriesname(96:113)); %HCP-D dtseries
    %uniqueID = char(dtseriesname(99:108)); %BCP8-11
    uniqueID = char(dtseriesname(100:109)); %BCP8-11
    
    disp(uniqueID)
    outputname = [outputdir filesep uniqueID];
    [~,whole_brain_number_of_nets,integration_zone_number_of_nets] = visualizedscalars(dtserieswithassignments{i},outputname,'number_of_networks',0);
    whole_brain_number_of_nets_all(i) = whole_brain_number_of_nets;
    integration_zone_number_of_nets_all(i) = integration_zone_number_of_nets;
    save([output_name '.mat'],'whole_brain_number_of_nets_all','integration_zone_number_of_nets_all');
end

