#!/usr/bin/env python

"""Conversion of .msh (tetrahedral) file to .ply surface files.
"""

import sys
import os
from os import chdir
import argparse
import numpy as np
from collections import OrderedDict
from pathlib import Path


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


def generate_ply_header(num_vertices, num_faces):
    """Helper function creating headers for .ply files.
    Args:
        num_vertices(int): The number of vertices.
        num_faces(int): The number of faces.
    
    Returns:
        ply_header(str): ply header for creation of ply file.
    """
    ply_header = """ply
format ascii 1.0
comment Created by meshio v5.0.0, 2021-08-11T19:09:09.536119
element vertex """+str(num_vertices)+"""
property double x
property double y
property double z
element face """+str(num_faces)+"""
property list uint8 int32 vertex_indices
end_header
"""
    return ply_header

def write_surfaces(filename, faces, vertices, dict_vertices):
    """Writes a .ply file based on vertices and faces and a
    dictionary containing information on how the vertices are connected.
    Performs a renumbering such that only essential vertices are output.

    Args: 
        filename(str): Name of the file.
        faces(np.ndarray): Array containing faces.
        vertices(np.ndarray): Array containing vertices.
        dict_vertices(dict): Dictionary for the given vertices. 
    
    Returns:
        None
    """
    num_faces = np.shape(faces)[0]
    num_vertices = np.shape(vertices)[0]

    # Print header to file
    outfile = open(filename, "w");
    ply_header = generate_ply_header(num_vertices, num_faces)
    outfile.write(ply_header)

    for i in range(0, num_vertices):
            outfile.write(str(vertices[i][0]) +" "+ str(vertices[i][1]) +" "+ str(vertices[i][2]) +"\n");

    # Remap vertices and faces
    for i in range(0, faces.shape[0]):
        vertex1 = faces[i][0]
        vertex2 = faces[i][1]
        vertex3 = faces[i][2]
        vertex1_new_id = dict_vertices[vertex1]
        vertex2_new_id = dict_vertices[vertex2]
        vertex3_new_id = dict_vertices[vertex3]
        outfile.write(f"3 {vertex1_new_id} {vertex2_new_id} {vertex3_new_id}\n")

    outfile.close()
    
def msh2arrays(name):
    
    """Returning array of nodes, vectors and elements from .msh file.
    
    Args:
        name(str): Name of the given file.

    Returns:
        node_data(np.ndarray): Array containing all the nodes/vertices.
        num_nodes(int): The number of nodes/vertices
        element_array_LV (np.ndarray): Array containing the faces on the left ventricle.
        element_array_RV (np.ndarray): Array containing the faces on the rightt ventricle.
        element_array_Epi (np.ndarray): Array containing the faces on the epicardium.
        element_array_MV (np.ndarray): Array containing the faces on the mitral valve.
        element_array_AV (np.ndarray): Array containing the faces on the aortic valve.
        element_array_PV (np.ndarray): Array containing the faces on the pulmonary valve.
        element_array_TV (np.ndarray): Array containing the faces on the tricuspid valve.
    """
    filename_in = name+'.msh'
    surfaces = {1:'LV', 2: 'RV', 3:'Epi', 4:'MV', 5:'AV', 6:'PV', 7:'TV'}
    element_array_LV = np.zeros((1,3), dtype=np.int32)
    element_array_RV = np.zeros((1,3), dtype=np.int32)
    element_array_Epi = np.zeros((1,3), dtype=np.int32)
    element_array_MV = np.zeros((1,3), dtype=np.int32)
    element_array_AV = np.zeros((1,3), dtype=np.int32)
    element_array_PV = np.zeros((1,3), dtype=np.int32)
    element_array_TV = np.zeros((1,3), dtype=np.int32)


    # read through the input file
    with open(filename_in, 'r') as f:
        lines = f.readlines()

    num_header_lines = 14 #manually adjust for now 
    line_pos = num_header_lines
    # read total number nodes 
    assert lines[line_pos].strip() == f"$Nodes"
    line_pos += 1
    # check that the number of nodes match our assumptions
    num_nodes = int(lines[line_pos].strip())
    line_pos += 1

    #save the coordinates to an array
    #remember meshalyzer/gmsh count from 1! but carp/paraview/python from 0!
    node_data = np.zeros((num_nodes, 3), dtype=np.float64)
    for i in range(num_nodes):
        line = lines[line_pos].strip()
        words = line.split()
        node_data[i, 0] = float(words[1])
        node_data[i, 1] = float(words[2])
        node_data[i, 2] = float(words[3])

        line_pos += 1
    #check that we are at the end of the Nodes
    assert lines[line_pos].strip() == f"$EndNodes"
    line_pos += 1 #to Elements
    assert lines[line_pos].strip() == f"$Elements"
    line_pos += 1 #to Elements
    # check element number total
    num_polygons_total = int(lines[line_pos].strip())
    polygon_data = np.zeros((num_polygons_total, 4), dtype=np.int32)
    line_pos += 1
    
    # filter element array for each surface using the surface dict
    while True:
        line = lines[line_pos].strip()
        words = line.split()
        if not((len(words) > 7) and (words[1] == '2')):
            line_pos+=1
        else:
            break
    assert (len(words) > 7) and (words[1] == '2') == True

    # polygon data start 
    # for each of the surfaces we create an array of elements 
    for key in surfaces:
        count_LV = 0
        count_RV = 0
        count_Epi = 0
        count_MV = 0
        count_AV = 0
        count_PV = 0
        count_TV = 0

        
        while True:
            line = lines[line_pos].strip()
            words = line.split()
            if ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_LV == 0) and (surfaces[key] == 'LV')):
                element_array_LV[0, 0] = int(words[5])-1
                element_array_LV[0, 1] = int(words[6])-1
                element_array_LV[0, 2] = int(words[7])-1
                line_pos += 1
                count_LV += 1
            
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_RV == 0) and (surfaces[key] == 'RV')):
                element_array_RV[0, 0] = int(words[5])-1
                element_array_RV[0, 1] = int(words[6])-1
                element_array_RV[0, 2] = int(words[7])-1
                line_pos += 1
                count_RV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_Epi == 0) and (surfaces[key] == 'Epi')):
                element_array_Epi[0, 0] = int(words[5])-1
                element_array_Epi[0, 1] = int(words[6])-1
                element_array_Epi[0, 2] = int(words[7])-1
                line_pos += 1
                count_Epi += 1
            
                    
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_MV == 0) and (surfaces[key] == 'MV')):
                element_array_MV[0, 0] = int(words[5])-1
                element_array_MV[0, 1] = int(words[6])-1
                element_array_MV[0, 2] = int(words[7])-1
                line_pos += 1
                count_MV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_AV == 0) and (surfaces[key] == 'AV')):
                element_array_AV[0, 0] = int(words[5])-1
                element_array_AV[0, 1] = int(words[6])-1
                element_array_AV[0, 2] = int(words[7])-1
                line_pos += 1
                count_AV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_PV == 0) and (surfaces[key] == 'PV')):
                element_array_PV[0, 0] = int(words[5])-1
                element_array_PV[0, 1] = int(words[6])-1
                element_array_PV[0, 2] = int(words[7])-1
                line_pos += 1
                count_PV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_TV == 0) and (surfaces[key] == 'TV')):
                element_array_TV[0, 0] = int(words[5])-1
                element_array_TV[0, 1] = int(words[6])-1
                element_array_TV[0, 2] = int(words[7])-1
                line_pos += 1
                count_TV += 1

            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_LV != 0) and (surfaces[key] == 'LV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_LV = np.vstack((element_array_LV,arr))
                line_pos += 1
                count_LV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_RV != 0) and (surfaces[key] == 'RV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_RV = np.vstack((element_array_RV,arr))
                line_pos += 1
                count_RV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_Epi != 0) and (surfaces[key] == 'Epi')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_Epi = np.vstack((element_array_Epi, arr))
                line_pos += 1
                count_Epi += 1
            
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_MV != 0) and (surfaces[key] == 'MV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_MV = np.vstack((element_array_MV, arr))
                line_pos += 1
                count_MV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_AV != 0) and (surfaces[key] == 'AV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_AV = np.vstack((element_array_AV, arr))
                line_pos += 1
                count_AV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_PV != 0) and (surfaces[key] == 'PV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_PV = np.vstack((element_array_PV, arr))
                line_pos += 1
                count_PV += 1
                
            elif ((len(words) > 7) and (words[1] == '2') and (words[3] == f'{key}') and (count_TV != 0) and (surfaces[key] == 'TV')):
                arr = np.zeros((1,3), dtype=np.int32)
                arr[0, 0] = int(words[5])-1
                arr[0, 1] = int(words[6])-1
                arr[0, 2] = int(words[7])-1
                element_array_TV = np.vstack((element_array_TV, arr))
                line_pos += 1
                count_TV += 1
                
            else:
                break
            
    # now we have all the node and element arrays
    return node_data,num_nodes, element_array_LV, element_array_RV, element_array_Epi, element_array_MV, element_array_AV, element_array_PV, element_array_TV
        
def filter_vertices_and_renumber(vertex_array_all, face_array_filtered):
    """Takes in arrays of vertices and faces and returns ordered dictionary
    of new vertex id and filtered vertices.
    
    Args:
        vertex_array_all(np.ndarray): Array containing all the vertices.
        face_array_filtered(np.ndarray): Array containing the filtered faces.

    Returns:
        filtered_nodes(np.ndarray): Array containing the filtered vertices.
        dicts(dict): Dictionary for the chosen vertices (new number to old one). 
    """
    
    values, counts = np.unique(face_array_filtered, return_counts=True)
    # filter vertices
    filtered_nodes = vertex_array_all[values]
    num_vertices = np.shape(filtered_nodes)[0]
    original_indices = values
    # create dictionary mapping new vertex numbers to old ones
    new_ids = list(range(0,num_vertices))
    original_ids_list = original_indices.tolist()
    d1 = zip(original_ids_list, new_ids)
    dicts = dict(d1)
    
    return filtered_nodes, dicts


if __name__ == '__main__':
    # get parent directory of CWD
    root = Path(".").resolve().parents[2]
    path2input = Path(root, 'data', 'input')
    path2out = Path(root, 'data', 'output')

    chdir(path2input)
    args = get_arguments(sys.argv)
    name = args.input_file

    # read in arrays
    node_data, num_nodes, element_array_LV, element_array_RV, element_array_Epi, element_array_MV, element_array_AV, element_array_PV, element_array_TV= msh2arrays(name)
    # reorder and write ply files 
    filtered_nodes_LV, dicts_LV = filter_vertices_and_renumber(node_data, element_array_LV)
    filtered_nodes_RV, dicts_RV = filter_vertices_and_renumber(node_data, element_array_RV)
    filtered_nodes_Epi, dicts_Epi = filter_vertices_and_renumber(node_data, element_array_Epi)
    filtered_nodes_MV, dicts_MV = filter_vertices_and_renumber(node_data, element_array_MV)
    filtered_nodes_AV, dicts_AV = filter_vertices_and_renumber(node_data, element_array_AV)
    filtered_nodes_PV, dicts_PV = filter_vertices_and_renumber(node_data, element_array_PV)
    filtered_nodes_TV, dicts_TV = filter_vertices_and_renumber(node_data, element_array_TV)
    
    new_out_directory = Path(path2out, str(name))
    chdir(new_out_directory)
    write_surfaces(f"{name}_endo_lv.ply", element_array_LV, filtered_nodes_LV, dicts_LV)
    write_surfaces(f"{name}_endo_rv.ply", element_array_RV, filtered_nodes_RV, dicts_RV)
    write_surfaces(f"{name}_epi.ply", element_array_Epi, filtered_nodes_Epi, dicts_Epi)
    write_surfaces(f"{name}_mv.ply", element_array_MV, filtered_nodes_MV, dicts_MV)
    write_surfaces(f"{name}_av.ply", element_array_AV, filtered_nodes_AV, dicts_AV)
    write_surfaces(f"{name}_pv.ply", element_array_PV, filtered_nodes_PV, dicts_PV)
    write_surfaces(f"{name}_tv.ply", element_array_TV, filtered_nodes_TV, dicts_TV)
    

    