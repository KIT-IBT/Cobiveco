#!/usr/bin/env python


"""Script containing several different methods used in extract_midmyocard.py
"""

import networkx as nx
import numpy as np
import re


def ply_reader(filename):
    """Reads the ply file to extract the vertices and the faces
    
    Args: 
        filename(str): The name of the given file.

    Returns:
        vertices(np.ndarray)
        faces(np.ndarray): Array containing the chosen faces.
        num_vertices(int): The number of vertices.
        num_faces(int): The number of faces.
    """
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
    return vertices, faces, num_vertices, num_faces

#unused here (?)
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

def write_surfaces(filename, face_tag_num, face_tag, faces, num_faces, vertices, num_vertices):
    """Outputs the .ply file with the extracted vertices and faces.
    Performs a renumbering such that only essential vertices are output.

    Args: 
        filename(str): The name of the chosen file.
        face_tag_num(int): Number of faces to output.
        face_tag(np.ndarray): Array containing the tags for the faces.
        faces(np.ndarray): Array containing the chosen faces.
        num_faces(int): The number of faces.
        vertices(np.ndarray): Array containing the chosen vertices.
        num_vertices(int): The number of vertices.
    
    Returns:
        None
    """
    # Count number of faces to output
    out_faces = 0
    out_vertices = 0
    vertex_tag = np.empty(num_vertices, dtype=int)
    vertex_tag.fill(-1)
    for i in range(0, num_faces):
        if face_tag[i] == face_tag_num:
            out_faces += 1
            for j in range(0, 3):
                vertex_tag[faces[i][j]] = 1

    # Renumber vertices
    renumber = 0
    for i in range(0, num_vertices):
        if vertex_tag[i] != -1:
            vertex_tag[i] = renumber
            renumber += 1
    out_vertices = renumber

    print("[info]", filename, "num_vertices: ", out_vertices)
    print("[info]", filename, "num_faces: ", out_faces)

    # Print header to file
    f2 = open(filename, "w");
    ply_header = generate_ply_header(out_vertices, out_faces)
    f2.write(ply_header)

    for i in range(0, num_vertices):
        if vertex_tag[i] != -1:
            f2.write(str(vertices[i][0]) +" "+ str(vertices[i][1]) +" "+ str(vertices[i][2]) +"\n");

    # Remap vertices and faces
    my_faces = np.zeros((faces.shape[0], faces.shape[1]), dtype=int)
    for i in range(0, faces.shape[0]):
        if face_tag[i] == face_tag_num:
            f2.write("3 ")
            for j in range(0, faces.shape[1]):
                f2.write(str(vertex_tag[faces[i][j]]) +" ");
            f2.write("\n")

    f2.close()

    
def find_face_id_numpy(vertex_pos, faces, num_faces):
    """Returns array of face ids for a given vertex

    Args: 
        vertex_pos(int): The position of the given vertex.
        faces(np.ndarray): Array containing the given faces.
        num_faces(int): The number of faces.

    Returns:
        matching_faces(np.ndarray): Array containing faces with a non-zero number of matches.
    """
    match_mask = faces == vertex_pos
    # reduce to 1D array with number of matching vertices for each face
    matching_vertex_per_face = np.sum(match_mask, axis=1)
    # get rows/faces with a non-zero number of matches
    matching_faces = np.flatnonzero(matching_vertex_per_face)
    return matching_faces