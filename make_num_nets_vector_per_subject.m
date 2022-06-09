function make_num_nets_vector_per_subject(dtseriesconc, outputdir)
%This function generates a vector od the number of networks for each
%subject based on a dtseries with overlapping network labels.

dtserieswithassignments = importdata(dtseriesconc);

for i = 1:length(dtserieswithassignments)
    disp(i)
    dtseriesname = dtserieswithassignments{i};
    %get output_names
    %HARDCODE WARNING
    %uniqueID = char(dtseriesname(88:110));
    %uniqueID = char(dtseriesname(101:123));
    uniqueID = char(dtseriesname(96:113)); %HCP-D dtseries
    disp(uniqueID)
    outputname = [outputdir filesep uniqueID];
    visualizedscalars(dtserieswithassignments{i},outputname,'number_of_networks',0)

end

