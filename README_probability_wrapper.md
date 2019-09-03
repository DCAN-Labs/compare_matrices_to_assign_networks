# Probability function wrapper

Python wrapper for `visualizedscalars` probability function script. This loads in dscalars to visualize them for all subjects. It also can calculate the probability of a network assignment from the provided list of subjects.

## Installation

Clone this repository and save it somewhere on the Linux/Unix system that you want to use it from.

## Dependencies

1. [Python version 3.5.2](https://www.python.org/downloads/release/python-352) or above
1. MathWorks [MATLAB Runtime Environment (MRE) version 9.1 (2016b)](https://www.mathworks.com/products/compiler/matlab-runtime.html)

## Usage

### Required Arguments

- `dscalarwithassignments` is one argument, the path to a `.dscalar.nii` file or to a `.conc` file which contains only a list of paths to `.dscalar.nii` files. If a .conc file is provided, it should have 1 path to a subject dscalar per row.
 
Here is an example of a basic call to this script:
```
python3 ./probability_function_wrapper.py example/file/path/dscalar_paths.conc
```

### Server-Dependent Arguments

- `--mre_dir` takes one argument, a path to the already-existing MATLAB Runtime Environment (MRE) directory. If this argument is excluded, then by default the wrapper will use hardcoded paths which are only valid on the RUSHMORE or Exacloud servers. If this script is run on a different server or locally, then this argument is required.

### Script Option Arguments

- `--DS_factor` takes one argument, a whole number (with only digits) to be used as the downsample factor. Reduce the 91282 vector by this factor to reduce the load on matlab visualization tools (e.g. `--DS_factor 2` will sample every other greyordinate, so visualization will have 45641 data points per subject).
- `--downsample_scalar` takes no arguments. If this flag is included, then the dscalar with be down sampled (sampled according to the `--DS_factor`). Otherwise, it will not downsample data for visualization.
- `--save_percentages` takes no arguments. Include this flag to save a file (dscalar for each network).

### Optional Path Arguments

- `--output` takes one argument, the path to a directory in which to place all of the files created by this script. If the directory does not exist, then it will be created at the given path. By default, this will be `./data`.
- `--tempfile` takes one argument, the path to a temporary file to write for saving.
- `--color_map` takes one argument, the path to file from which to load RGB values for `colormap` to match network colors.

### Advanced Usage

Here is an example of a more complex call to this script, to create visualizations for all dscalars listed in `data_dscalar_paths.conc` and downsample so that only every 50th greyordinate is sampled. Once the script is finished, it will save the `dscalar` output files to the `./visualizations` directory.
```
python3 ./probability_function_wrapper.py ./example/file/path/data_dscalar_paths.conc --downsample_scalar --DS_factor 50 --save_percentages --output ./visualizations
```

For more usage information, call `python3 probability_function_wrapper.py --help`.

## Required Files

- `./probability_function_wrapper.py`
- `./support_files/PowerColorMap.mat`
- `./support_files/run_visualizedscalars.sh`
- `./support_files/visualizedscalars`

## Updates

2019-08-30: Documentation initially written by Greg Conan.