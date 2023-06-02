#!/usr/bin/env python3

"""Filtering triangles, that are wrongly assumed
part of the bridge/boundary.
Written by: Aadarsh Bussooa, Lisa Pankewitz
"""

import networkx as nx
import numpy as np
import re
from extractfunc_bridge import *
import argparse
import sys


def get_arguments(input_args):
    """Returns the name of the file on which the script is run.

    Args:
        input_args(list): List containing the name of the script and the file on which it runs.

    Returns:
        parser.parse_args() (Namespace): The name of the file.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("initial_triangle",type=int, help='Enter input name of .pts .elem .surf files without extension')
    parser.add_argument("input_file_name", help='Enter output name without extension')
    parser.add_argument("output_file_name", help='Enter output name without extension')
    return parser.parse_args()


if __name__ == "__main__":
    
    args = get_arguments(sys.argv)
    vertex_initial = args.initial_triangle
    filename = args.input_file_name
    output = args.output_file_name
    
    # set in parameter file
    # Heuristics (should usually not be changed)
    base_depth_limit = 1 # identifies only surronding nodes around seed
    base_iterations = 20 # iterate multiple times to identify all connected nodes
    surface_threshold = 0.25 # 0.25 rad; 11 deg
    sensitive_threshold = 0.1 # 0.1 rad; 5.7 deg
    surface_depth_limit = 550 # used for epi and endo
    
    # Read file 
    f = open(filename, "r")
    lines = f.readlines()

    # Read ply header
    end = 10
    for i in range(0, end):
        keyword = re.search("^element vertex (\d*)", lines[i])
        if(keyword != None):
            num_vertices = int(keyword.group(1))
        keyword = re.search("^element face (\d*)", lines[i])
        if(keyword != None):
            num_faces = int(keyword.group(1))

    print("[info] num_vertices: ", num_vertices)
    print("[info] num_faces: ", num_faces)

    # Read vertices
    start = 10
    vertices = np.zeros((num_vertices, 3), dtype=float)
    for i in range(start, start+num_vertices):
        data = re.search("^([e0-9.-]*) ([e0-9.-]*) ([e0-9.-]*)", lines[i])
        vertices[i-start] = [float(data.group(1)), float(data.group(2)), float(data.group(3))]

    # Read faces
    start = 10 + num_vertices
    faces = np.zeros((num_faces, 3), dtype=int)
    for i in range(start, start+num_faces):
        data = re.search("^3 ([0-9]*) ([0-9]*) ([0-9]*)", lines[i])
        faces[i-start] = [int(data.group(1)), int(data.group(2)), int(data.group(3))]

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
    id_1 = find_face_id(vertex_initial, faces, num_faces)


    # Apply BFS algorithm to find connected vertices of the base
    # Tag definition: 0 untouched; >0 BFS applied; -1 disabled
    # Seed point initially tagged with 1
    print("[info] extracting base")
    face_tag = np.zeros(num_faces, dtype=int)
    face_tag_sum = face_tag
    face_tag[id_1[0]] = 1
    face_tag, tag = apply_bfs_edge_detect(face_tag, faces, num_faces, g, base_depth_limit, vertices, base_iterations, surface_threshold)
    write_surfaces(f"{output}.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    face_tag_sum += face_tag