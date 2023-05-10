#!/usr/bin/env python3

"""Automated surface file generation from cardiac meshes
Creates base, epi and endo surfaces from ply mesh
Contains commented function to extract midmyocardium
Written by: Aadarsh Bussooa, Lisa Pankewitz
Affiliation: Simula Research Laboratory
"""

import networkx as nx
import numpy as np
import re
from extractfunc import *


import argparse
import sys


def get_arguments(input_args):
    """Argument parser. Returns list of arguments passed as commandline arguments.

    Args:
        input_args(list): List of input arguments.

    Returns:
        parser.parse_args() (Namespace): The name of the file.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_epicardium",help='Enter input name of _epi.ply file with extension')
    return parser.parse_args()


if __name__ == "__main__":
    
    args = get_arguments(sys.argv)

    # Parameters
    filename = args.input_epicardium
    # file base
    output = filename[:-8]
    # inputfile names of base 
    base1 = output+"_mv.ply"
    base2 = output+"_av.ply"
    base3 = output+"_tv.ply"
    base4 = output+"_pv.ply"

    # Heuristics (should usually not be changed)
    base_depth_limit = 1 # identifies only surronding nodes around seed
    base_iterations = 10 # iterate multiple times to identify all connected nodes
    surface_threshold = 0.25 # 0.25 rad; 11 deg
    sensitive_threshold = 0.1 # 0.1 rad; 5.7 deg
    surface_depth_limit = 550 # used for epi and endo
    axis_up = [1, 0, 0]
    inclination_threshold = 0.55

    vertices, faces, num_vertices, num_faces = ply_reader(filename)
    vertices_base1, faces_base1, num_vertices_base1, num_faces_base1 = ply_reader(base1)
    vertices_base2, faces_base2, num_vertices_base2, num_faces_base2 = ply_reader(base2)
    vertices_base3, faces_base3, num_vertices_base3, num_faces_base3 = ply_reader(base3)
    vertices_base4, faces_base4, num_vertices_base4, num_faces_base4 = ply_reader(base4)

    # calculate starting vertex
    centroid = np.mean(np.concatenate((vertices_base1,vertices_base2,vertices_base3,vertices_base4,vertices_base4),axis=0),0)

    # Construct graph from face data
    print("[info] constructing graph")
    g = nx.Graph()
    for i in range(0, num_faces):
        g.add_edge(int(faces[i][0]), int(faces[i][1]))
        g.add_edge(int(faces[i][1]), int(faces[i][2]))
        g.add_edge(int(faces[i][2]), int(faces[i][0]))
        g.nodes[int(faces[i][0])][i] = "key"
        g.nodes[int(faces[i][1])][i] = "key"
        g.nodes[int(faces[i][2])][i] = "key"

    # Find face_id from vertex

    face_tag = np.zeros(num_faces, dtype=int)

    # calculate initial centroid from the center between all the valves!
    distances = np.linalg.norm(vertices - centroid, axis=1)
    indices = []
    vertices_to_be_excluded = np.where(distances < 26)
    for vertex in range(0, np.shape(vertices_to_be_excluded)[1]):
        index = find_face_id_numpy(vertices_to_be_excluded[0][vertex], faces, num_faces)
        indices.append(index)
    flattened_indices = (list(np.concatenate(indices).flat))
    unique_indices = np.unique(flattened_indices)
    # find face_id
    face_tag = np.zeros(num_faces, dtype=int)

    face_tag[unique_indices] = 1

    write_surfaces(f"{output}_epi_base.ply", 1, face_tag, faces, num_faces, vertices, num_vertices)
    write_surfaces(f"{output}_epi_no_base.ply", 0, face_tag, faces, num_faces, vertices, num_vertices)
