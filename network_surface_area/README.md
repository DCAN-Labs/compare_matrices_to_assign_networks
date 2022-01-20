# Network Surface Area Documentation 
### Creating a dense surface area map from a given neural network 

The human brain can tell us many things about its functionality, but a lot of it is obscured if analyzed incorrectly, or specifically, if brain volume is used to analyze network connectivity across brain regions. This lack of accuracy is due to the fact that the volume does not account for the many folds in the brain, misproportionatly shortening their distance when viewed in Euclidean distance. Because of this, surface area is preferred. This code creates a surface area map of 14 neural networks (Gordon et al. 2017) using left and right mid-thickness files of the brain. 

This package will create a dense scalar map in CIFTI format (dscalar.nii) [(Glasser et al. 2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/) showing surface area from given left and right midthickness files. 

## Where is the code?

This procedure is part of a package called Compare Matrices to Assign Networks.

The surface area package can be found at: https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks/-/blob/master/network_surface_area_from_network_file.m

***

## Background

The input data are assumed to be processed with the FreeSurfer processing pipeline. The code assumes that particpants have 91282 greyordinates cortex + subcortical structures or are cortex only 549412.

The code works by taking a dense scalar matrix (dscalar.nii) and uses the left and right mid-thickness files (.gii) provided to create a map of surface area map based on 14 pre-defined neural networks ([Gordon et al. 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X?via%3Dihub)). 

[!Placeholder image: MSC average]

This documentation gives an example using data from the Midnight Scan Club (MSC) which can be found [here](https://github.com/MidnightScanClub/MSCcodebase/tree/master/Utilities/Conte69_atlas-v2.LR.32k_fs_LR.wb). 

## Create files ahead of time

### Required:
**Lmidthicknessfile -** a left mid-thickness file of the brain (.surf.gii)

**Rmidthicknessfile -** a right mid-thickness file (.surf.gii)

**dscalarwithassignnments -** vector of assignments for neural network numbers 1-14


## Example Call

The most current iteration of the code is *network_surface_area_from_network_file.m*

[!Placeholder Image (2)]


`network_surface_area_from_network_file(Lmidthicknessfile,Rmidthicknessfile,dscalarwithassignments,outputname)
`

- **Lmidthicknessfile =** a left mid-thickness file of the brain
- **Rmidthicknessfile =** a right mid-thickness file of the brain

- **dscalarwithassignnments =** vector of assignments for neural network numbers 1-14
- **outputname =** name of resulting file 

### Output:
- dscalar of surface area of the left and right brain, overlayed with neural networks data (.nii)


