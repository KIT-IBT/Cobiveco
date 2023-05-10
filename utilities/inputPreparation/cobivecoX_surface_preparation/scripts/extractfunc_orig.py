#!/usr/bin/env python


"""Script containing several different methods used in extract.py
"""

import networkx as nx
import numpy as np
import re


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
    f2 = open(filename, "a");
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

#maybe say what bfs stands for?
def apply_bfs_algorithm(face_tag, faces, num_faces, g, bfs_depth_limit):
    
    """Applies BFS with higher depth to quickly identify connected nodes
    Used to identify surfaces after removing the base

    Args:
        face_tag(np.ndarray): Array containing the initial tags for the faces.
        faces(np.ndarray): Array containing the chosen faces.
        num_faces(int): The number of faces.
        g(Graph): The given graph.
        bfs_depth_limit(int): The upper bound on the depth.

    Returns:
        face_tag(np.ndarray): Array containing the new tags for the faces.
        tag(int): Tag for the faces.
    """
    tag = 1
    bfs_results = np.empty((0), dtype=int)
    bfs_results_keys = np.empty(0, dtype=int)
    for i in range(0, num_faces):
        if face_tag[i] == tag:
            face_tag[i] = tag+1
            try:
                bfs_results = np.array(nx.bfs_tree(g, source=faces[i][0], depth_limit=bfs_depth_limit).nodes())
            except Exception:
                pass
            for j in bfs_results:
                bfs_results_keys = list(g.nodes[j].keys())
                for k in bfs_results_keys:
                    if(face_tag[k] != -1):
                        face_tag[k] = tag+1
    return face_tag, tag+1


def find_face_id(vertex_pos, faces, num_faces):
    """Returns array of face ids for a given vertex.

    Args: 
        vertex_pos(int): Id of the starting vertex.
        faces(np.ndarray): Array containing the chosen faces.
        num_faces(int): The number of faces.

    Returns:
        search_res(np.ndarray): Array with all the faces' ids.
    """
    search_res = np.array([], dtype=int)
    for i in range(0, num_faces):
        if any(faces[i] == vertex_pos):
            search_res = np.append(search_res, i)
    return search_res

