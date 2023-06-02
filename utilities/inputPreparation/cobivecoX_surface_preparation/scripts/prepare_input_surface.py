#!/usr/bin/env python

"""Batch processing of input file creation.
"""

import subprocess
from pathlib import Path  
from os import chdir
import shutil
import pandas as pd
import numpy as np 
import csv
import configparser


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

def createconfigtemplate(file):
    """Helper function to create the template for the configuration file.

    Args:
        file(str): Name of file.

    Returns:
        configfile(TextIOWrapper): Template to the configuration file needed to extract the surfaces.
    """
    config = configparser.ConfigParser()
    config['start_indices'] = {'mv': '',
                               'av': '',
                               'pv': '',
                               'tv': '',
                               'epi': '',
                               'rv': '',
                               'lv': ''}
    with open(file+'.config', 'w') as configfile:
        config.write(configfile)
    return configfile

if __name__ == "__main__":
    # get parent directory of CWD
    root = Path(".").resolve().parent
    path2input = Path(root, 'data', 'input')
    path2out = Path(root, 'data', 'output')
    # File that contains a list of files names 
    file = Path(path2input, 'example.xlsx')

    # get ids from all files/dir/models you want to look into
    id_array = readatlas(file)
    
    # collect names of all meshes that successfullly created in the list
    ids_created = []
    ids_failed = []
    
    # Make a data output directory for each mesh ID
    # make a gmsh directory
    for i in range(len(id_array)):
        # go to input directory and check that the right files are there
        # Make a directory with ID name for each mesh where we collect
        # all the output files
        dir_name = id_array[i]
        chdir(path2input)
        # check that files exist
        meshname = str(dir_name)+".msh"
        meshname_vol = str(dir_name)+"_vol.msh"
        # create path 
        path_to_meshfile = Path(path2input, str(meshname))
        path_to_meshfile_vol = Path(path2input, str(meshname_vol))
        
        if path_to_meshfile.stat().st_size < 500 or (path_to_meshfile.exists() is not True) or path_to_meshfile_vol.stat().st_size < 500 or (path_to_meshfile_vol.exists() is not True):
            print(f"[info] Something went wrong with {dir_name}. Either failed volume mesh, mixed mesh or both. Please check. Adding to the failed list.")
            # add to failed list
            ids_failed.append(dir_name)
        elif path_to_meshfile.stat().st_size > 5000 and (path_to_meshfile.exists() is True) and path_to_meshfile_vol.stat().st_size > 5000 and (path_to_meshfile_vol.exists() is True):
            ids_created.append(dir_name)
            # mkdir for the outputfiles
            new_out_directory = Path(path2out, str(dir_name))
            # change dir to make the files in the right spot
            if not Path(new_out_directory).exists():
                new_out_directory.mkdir(parents=True, exist_ok=True)
            

            with open(meshname, 'r') as f:
                lines = f.readlines()
                if ('$PhysicalNames\n' and '8\n') in lines: 
                    # if the gmsh file is labeled, extract surfaces using those labels
                    chdir(new_out_directory)
                    list_commands = subprocess.call([f"python ../../../scripts/msh2surfacesply.py {dir_name}"], shell=True)
                else:
                    # if the gmsh file is not labeled, first convert it to ply
                    chdir(new_out_directory)
                    list_commands = subprocess.call([f"python ../../../scripts/msh2ply.py {dir_name}"], shell=True)

                    path_to_config = Path(dir_name+'.config')
                    # create a config template if the config file is not provided
                    if not path_to_config.exists():
                        createconfigtemplate(dir_name)

                    with open(dir_name+'.config', 'r') as f:
                        # read the config file
                        lines = f.readlines()

                        line_pos = 0
                        # does not go on if [start_indices] is not in config file (or in wrong position)
                        assert lines[line_pos].strip() == f"[start_indices]"
                        line_pos += 1
                        # print an error message if the config file does not contain the required ids
                        for i in range(line_pos,7):
                            line = lines[i].strip()
                            words = line.split()
                            assert len(words) == 3, "The configuration file does not contain the required start indices, please modify it according to the README file"
                            i += 1

                        # extract the surfaces from the ply using the given start indices
                        list_commands = subprocess.call([f"python ../../../scripts/extract.py {dir_name}.ply"], shell=True)

            # Extract midmyocard
            list_commands = subprocess.call([f"python ../../../scripts/extract_midmyocard.py {dir_name}_epi.ply"], shell=True)

            ids_created.append(dir_name)
        else:
            print("Something went wrong here.")

    # Copy files to directory in examples
    list_commands = subprocess.call([f"python ../../../scripts/cp2cobiveco.py"], shell=True)

    # export csv of failed ids
    chdir(path2out)
    ids_created_a = np.asarray(ids_created)
    ids_failed_a = np.asarray(ids_failed)

    np.savetxt("ids_created_inputpreparation.csv",ids_created_a, delimiter=',', fmt="%s")
    np.savetxt("ids_failed_inputpreparation.csv",ids_failed_a, delimiter=',', fmt="%s")