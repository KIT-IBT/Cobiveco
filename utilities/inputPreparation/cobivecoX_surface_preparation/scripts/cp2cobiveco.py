#!/usr/bin/env python

"""Copy the needed cobivecox input files to the right directory.
"""

import subprocess
from pathlib import Path  
from os import chdir
import shutil
import pandas as pd
import numpy as np 
import csv



def readatlas(file):
    """Returns an array with the ids of all the files you want to look into.
    
    Args:
        file(Path): Excel file containing files names.
    
    Returns:
        array(np.ndarray, str): Array with strings of filenames.
    """
    df = pd.read_excel(file)
    array = df['Instance ID'].to_numpy(dtype=str)
    # iterate over array[i] for ID
    return array

if __name__ == "__main__":
    root = Path(".").resolve().parents[2]
    path2input = Path(root, 'data', 'input')
    path2out = Path(root, 'data', 'output')
    file = Path(path2input, 'example.xlsx')

    
    
    # get ids from all files/dir/models you want to look into
    id_array = readatlas(file)
    
    ids_created = []
    ids_failed = []
    
    # Make a data output directory for each mesh ID
    # make a gmsh directory
    for i in range(len(id_array)):
        # got to gmsh directory and check that msh files were created
        
        # Make a directory with ID name for each mesh where we collect
        # all the output files
        dir_name = id_array[i]
        new_out_directory = Path(path2out, str(dir_name))

        path2cobiveco = Path(".").resolve().parents[5]
        copy_directory = Path(path2cobiveco, 'examples', str(dir_name))
        if not Path(copy_directory).exists():
            copy_directory.mkdir(parents=True, exist_ok=True)
        # Make a directory for the gmsh in- and output files
        new_gmsh_directory = Path(new_out_directory, 'gmsh')

        # check that files exist
        meshname_vol = str(dir_name)+".msh"
        meshname = str(dir_name)+"_vol.msh"
        # create path 
        path_to_meshfile = Path(new_gmsh_directory, str(meshname))
        path_to_meshfile_vol = Path(new_gmsh_directory, str(meshname_vol))
        # mkdir for cobiveco input file preparation
        inputfile_prep_path = Path(new_gmsh_directory, 'inputfile_prep')
        if not Path(inputfile_prep_path).exists():
            inputfile_prep_path.mkdir(parents=True, exist_ok=True)
        # remesh files
        # change dir to make the files in the right spot
        chdir(new_out_directory)
        
        list_commands = subprocess.call([f"cp {dir_name}.vtu  ../../../../../../examples/{dir_name}/"], shell=True)
        list_commands = subprocess.call([f"cp {dir_name}_epi_*.ply ../../../../../../examples/{dir_name}/"], shell=True)
        list_commands = subprocess.call([f"cp {dir_name}_endo_*.ply ../../../../../../examples/{dir_name}/"], shell=True)
        list_commands = subprocess.call([f"cp {dir_name}_*v.ply ../../../../../../examples/{dir_name}/"], shell=True)

    # export csv of failed ids
    # change directory to 
    chdir(path2out)
    ids_created_a = np.asarray(ids_created)
    ids_failed_a = np.asarray(ids_failed)
    np.savetxt("ids_created_cp.csv",ids_created_a, delimiter=',', fmt="%s")
    np.savetxt("ids_failed_cp.csv",ids_failed_a, delimiter=',', fmt="%s")