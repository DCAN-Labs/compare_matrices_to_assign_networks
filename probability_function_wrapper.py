#! /usr/bin/env python3

"""
Probability Function Wrapper
Greg Conan: conan@ohsu.edu
Created 2019-07-23
Last Updated 2019-08-30
"""

##################################
#
# Wrapper for probability function script which loads in a .conc of dscalars
# to visualize them for all subjects.
#
##################################

import argparse
import os
from socket import gethostname
import subprocess

# Constants

# Default path to temporary file, color map file, and output directory
DEFAULT_COLORMAP = "./support_files/PowerColorMap.mat"
DEFAULT_OUTPUT = "./data"
DEFAULT_TEMPFILE = "./support_files/test_clean_91282.dscalar.nii"

# MATLAB Runtime Environment (MRE) directories which depend on the host server
MRE_EXACLOUD = ("/home/exacloud/lustre1/fnl_lab/code/external/utilities/"
                "Matlab2016bRuntime/v91")
MRE_RUSHMORE = ("/mnt/max/shared/code/external/utilities/Matlab2016bRuntime/"
                "v91")

# BASH script to run compiled MATLAB visualizedscalars script
SCRIPT_VISUALIZEDSCALARS = "./support_files/run_visualizedscalars.sh"


def main():
    cli_args = get_cli_args()

    print("\nProbability Function Wrapper started running at: ")
    subprocess.check_call("date")

    visualizedscalars(cli_args)

    print("\nProbability Function Wrapper finished at: ")
    subprocess.check_call("date")


def get_cli_args():
    """
    Get and validate all args from command line using argparse.
    :return: Namespace with all validated inputted command line args
    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        description=("Wrapper for probability function script. This loads in "
                     "a conc of dcalars to visualize them for all subjects. "
                     "It also can calculate the probability of a network "
                     "assignment from the provided list of subjects.")
    )

    # Required: Path to dscalar input file
    parser.add_argument(
        "dscalarwithassignments",
        help=("Path to a .dscalar.nii file, or to a .conc file which contains "
              "only a list of paths to .dscalar.nii files. If a .conc file is "
              "provided, it should have 1 path to a subject dscalar per row.")
    )

    # Optional: Output folder
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=("Directory to place all of the files created by this script. If "
              "the directory does not exist, then it will be created at the "
              "given path. By default, this will be " + DEFAULT_OUTPUT)
    )

    # Optional parameters originally from MATLAB script
    parser.add_argument(
        "-D",
        "--DS_factor",
        type=valid_whole_number_as_string,
        help=("Downsample factor. Reduce the 91282 vector by this factor to "
              "reduce the load on matlab visualization tools (e.g. DS = 2 "
              "will sample every other greyordinate, so visualization will "
              "have 45641 data points per subject). This must be a whole "
              "number (with only digits).")
    )
    parser.add_argument(
        "-d",
        "--downsample_scalar",
        default="0",
        action="store_const",
        const="1",
        help=("If this flag is included, then the dscalar with be down "
              "sampled (sampled according to the DS_factor). Otherwise, it "
              "will not downsample data for visualization.")
    )
    parser.add_argument(
        "-s",
        "--save_percentages",
        default="0",
        action="store_const",
        const="1",
        help="Include this flag to save a file (dscalar for each network)."
    )

    # Optional: Get path to MRE / MCR root to run compiled MATLAB script
    parser.add_argument(
        "-m",
        "--mre_dir",
        type=str,
        help=("Path to the directory containing the MATLAB Runtime "
              "Environment (MRE) version 9.1 or newer. This is used to run "
              "a compiled MATLAB script. This argument must be a valid path "
              "to an existing folder.")
    )

    # Optional: open a temp file to write for saving
    parser.add_argument(
        "-t",
        "--tempfile",
        default=DEFAULT_TEMPFILE,
        help=("Path to temporary file to write for saving. Default file path "
              "is " + DEFAULT_TEMPFILE)
    )

    # Optional: Get color map file path
    parser.add_argument(
        "-c",
        "--color_map",
        default=DEFAULT_COLORMAP,
        help=("Path to file from which to load RGB values for colormap to "
              "match network colors. Default path is " + DEFAULT_COLORMAP)
    )

    # Parse, validate, and return all CLI args
    return validate_cli_args(parser.parse_args(), parser)


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    # Validate DS factor
    if int(cli_args.DS_factor) >= 91282:
        parser.error("Downsample factor must be less than 91282.")

    # Validate file/dir path CLI args
    validate_readable_file(cli_args.dscalarwithassignments, parser)
    validate_readable_file(cli_args.color_map, parser)
    validate_readable_file(cli_args.tempfile, parser)
    validate_and_get_output_dir(cli_args.output, parser)

    # If no MRE dir was provided, use a default depending on the host server
    if not cli_args.mre_dir:
        host = gethostname()
        if host == "rushmore":
            cli_args.mre_dir = MRE_RUSHMORE
        elif "exa" in host:
            cli_args.mre_dir = MRE_EXACLOUD
        else:
            parser.error("Please enter a path to a valid MATLAB Runtime "
                         "Environment directory in the --mre_dir argument.")
        if not os.access(cli_args.mre_dir, os.R_OK):
            parser.error("Cannot read MATLAB Runtime Environment directory at "
                         + cli_args.mre_dir)

    return cli_args


def validate_readable_file(path, parser):
    """
    Throw exception unless parameter is a valid readable filename string. This
    is used instead of argparse.FileType("r") because the latter leaves an open
    file handle, which has caused problems.
    :param path: Parameter to check if it represents a valid filename
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: A valid filename as a string
    """
    try:
        assert os.access(path, os.R_OK)
        return os.path.abspath(path)
    except (AssertionError, OSError, TypeError):
        parser.error("Cannot read file at " + path)


def valid_whole_number_as_string(str_to_check):
    """
    Throws ValueError unless str_to_check is a string representing a positive
    integer.
    :param str_to_check: string to validate
    :return: string which is either "none" or represents a valid float value
    """
    assert str_to_check.isdigit()
    return str_to_check


def validate_and_get_output_dir(output, parser):
    """
    Validate that output directory exists and is writeable, or create one, and
    then return valid path to that output directory without trailing slash
    :param output: File path of output directory to validate
    :param parser: argparse ArgumentParser to raise error if path is invalid
    :return: Valid path to output directory
    """
    try:
        output = os.path.abspath(output)
        if not os.path.isdir(output):
            os.makedirs(output)
        assert os.access(output, os.W_OK)
        return output[:-1] if output[-1] == "/" else output
    except (OSError, TypeError):
        parser.error("Cannot make output folder at " + output)
    except AssertionError:
        parser.error("Cannot write to output folder at " + output)


def visualizedscalars(cli_args):
    """
    Run compiled MATLAB script to visualize dscalars
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    subprocess.check_call((
        SCRIPT_VISUALIZEDSCALARS,
        cli_args.mre_dir,
        cli_args.dscalarwithassignments,
        cli_args.downsample_scalar,
        cli_args.DS_factor,
        cli_args.save_percentages,
        cli_args.tempfile,
        cli_args.color_map,
        cli_args.output
    ))


if __name__ == '__main__':
    main()
