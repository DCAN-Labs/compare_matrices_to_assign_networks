# How to obtain Border Perimeters
### Calculating lengths of neural networks for further analysis

## Border Perimeter Importance
To calculate inter-greylength in the brain, surface files and their respective data, including border perimeters of neural networks, are necessary. 

## Code Location
This code is part of a package called Compare Matrices to Assign Networks, found [here](https://gitlab.com/Fair_lab/compare_matrices_to_assign_networks).

## Step-by-Step Tutorial
This tutorial uses data from the Midnight Scan Club (MSC) to demonstrate the code's function. It can be found [here](https://github.com/MidnightScanClub/MSCcodebase). Specifically, we are using data from the folder [Conte69_atlas-v2.LR.32k_fs_LR.wb](https://github.com/MidnightScanClub/MSCcodebase/tree/master/Utilities/Conte69_atlas-v2.LR.32k_fs_LR.wb), under the Utilities file, within the application workbench. 

### Background
This code requires six inputs described below. A majority of these input will require *.conc files, which are files that translate information with *Concurnas* Source Code into list format. For example, a dlabel file, a dense map of numerical values, will be a list of paths *to* the desired dlabel files in *.conc format. 

1. **`dlabel_conc`** = name of *.conc file containing a list of paths to the dlabel files of the data you want to analyze 
    * *converting a dscalar file to a dlabel file:* use "wb_command-cifti-label-import mydscalar.nii mylabeltable.txt mydlabel.nii"

2. **`L_surf_conc`** = name of *.conc file containing a list of paths to the left midthickness surface area data

3. **`R_surf_conc`** = name of *.conc file containing a list of paths to the right midthickness surface area data
    * *important note:* make sure surface area data values are consistent with the dlabel file values  

4. **`remove_border_files`** = 0 or 1
    * to calculate border perimeters, this code creates files of borders based on the [Glasser et al. (2013)](https://pubmed.ncbi.nlm.nih.gov/23668970/) parcellations
        * set to 0 to *keep* these files after the code has run 
        * set to 1 to *remove* these files after the code has run

5. **`get_compactness`** = 0 or 1
    * this part(/input??) will calculate network compactness using the Polsby-Popper test and provide them as an output file 
        * set to 0 for only the border perimeters
        * set to 1 for the border perimeters and the compactness values

6. **`varargin`** = name of *.mat file or *.conc file containing a list of paths to surface area data
    * *important note:* if it is a *.conc file, make sure it is a list of *.label files consistent with the values of the above *.conc files and in the same order



### Example Call
`getborderperimeters(dlabel_conc,L_surf_conc,R_surf_conc, remove_border_files,get_compactness,varargin)`

### Step ?: Running locally vs. separately(?) 

### Step 0: Validating input file existence
* Lines 64-114  [maybe replace with pics for visual? still many pictures...] of this code double check to make sure the input files and links provided exist and can be used to obtain border perimeters
    * Order of validation is (1) *.dtseries; (2) *.dscalar; and (3) *.motion(?) 
    * If the input *.conc files' values are not consistent with each other, an error message will be displayed stating that resulting subject data will be out of order (Notes say code will add empty cells to aid correlation?)


### Step 1: Checking input file type is dlabel
* Similar to Step 0, Step 1 uses lines 132-136 to make sure the input file can be used for this code, in this case, checking to make sure the file type has dlabel data

### Step 2: Seperating dlabels into two *.giis
* Lines 139-141 will seperate the dlabel data files into two GIFTI-formatted files 

### Step 3: Converting labels to borders
* Lines 144-149 will convert the dlabel data, now in GIFTI format, into border data ![Labeled border](images/perimeter_illus1.jpg?raw=true "Dlabel files are converted to border files")

### Step 4: Calculating border length
* When the data is finally border data, lines 153-239 will calculate the length of the borders
    * This will occur with the left hemisphere first and then the right
    * For each hemisphere calculation, an empty table will be created to fill respectively with data corresponding to the original *.conc lists' data 
    * The results of this will be separate left and right tables and a scatter plot of the two plotted against the other

#### Step 4a: (Optional) Removing border files 
* If requested (`remove_border_files=1`), lines 228-239 will delete the intermediate border files created to calculate the perimeters 


## Outputs
There will be 2-5 output files depending on what was entered as the input. 

* 2 tables of border length values of data provided (*.conc(??))
    * One left, one right
* *If requested:* 2 border files (*.border)
    * One left, one right
* *If requested:* 1 network compactness values file (*.conc)
