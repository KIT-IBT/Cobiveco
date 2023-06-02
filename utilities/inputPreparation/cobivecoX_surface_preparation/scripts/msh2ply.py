#!/usr/bin/env python

"""Batch conversion of .msh files to .vtu files and .ply file
"""

import sys
import subprocess
from pathlib import Path  # instead of os module
from os import chdir
import shutil
import pandas as pd
import numpy as np 
import csv
import argparse


def get_arguments(input_args):
    """Argument parser. Returns list of arguments passed as commandline arguments.

    Args:
        input_args(list): List of input arguments.

    Returns:
        parser.parse_args() (Namespace): The name of the file.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_file", help='Enter input name of .pts .elem .surf files without extension')
    return parser.parse_args()

def readatlas(file):
    """Returns an array of strings from an excel file with file names/ids.
    
    Args:
        file(Path): Path to excel file containing file names.
    
    Returns:
        array(np.ndarray, str): Array containing file names.
    """
    df = pd.read_excel(file)
    array = df['Instance ID'].to_numpy(dtype=str)
    # iterate over array[i] for ID
    return array

if __name__ == '__main__':
    # get parent directory of CWD
    root = Path(".").resolve().parents[2]
    path2input = Path(root, 'data', 'input')
    path2out = Path(root, 'data', 'output')

    chdir(path2input)
    args = get_arguments(sys.argv)
    dir_name = args.input_file
    # check that files exist
    meshname = str(dir_name)+".msh"
    meshname_vol = str(dir_name)+"_vol.msh"
    # create path 
    path_to_meshfile = Path(path2input, str(meshname))
    path_to_meshfile_vol = Path(path2input, str(meshname_vol))
        
    if path_to_meshfile.stat().st_size < 500 or (path_to_meshfile.exists() is not True) or path_to_meshfile_vol.stat().st_size < 500 or (path_to_meshfile_vol.exists() is not True):
        print(f"[info] Something went wrong with {dir_name}. Either failed volume mesh, mixed mesh or both. Please check. Adding to the failed list.")
    elif path_to_meshfile.stat().st_size > 5000 and (path_to_meshfile.exists() is True) and path_to_meshfile_vol.stat().st_size > 5000 and (path_to_meshfile_vol.exists() is True):
        # mkdir for the outputfiles
        new_out_directory = Path(path2out, str(dir_name))
        # change dir to make the files in the right spot
        if not Path(new_out_directory).exists():
            new_out_directory.mkdir(parents=True, exist_ok=True)
        # change directory to make the files in the right location
        # create .ply file using meshio
        list_commands = subprocess.call([f"meshio convert {dir_name}.msh ../output/{dir_name}/{dir_name}.ply --ascii"], shell=True)
        # create .vtu file using meshio
        list_commands = subprocess.call([f"meshio convert {dir_name}_vol.msh ../output/{dir_name}/{dir_name}.vtu"], shell=True)
    else:
        print("Something went wrong here.")
    