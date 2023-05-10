#!/usr/bin/env python

'''Automated surface file generation from cardiac meshes
Creates base, epi and endo surfaces from ply mesh when you do not have them labeled
Written by: Aadarsh Bussooa (2021), Lisa Pankewitz (2022)
'''

import networkx as nx
import numpy as np
import re
from extractfunc_orig import *
import argparse
import sys
import configparser
import queue


def get_arguments(input_args):
    """Returns the name of the file on which the script is run.

    Args:
        input_args(list): List containing the name of the script and the file on which it runs.

    Returns:
        parser.parse_args() (Namespace): The name of the file.
    """

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_filename", help='Enter input name of .pts .elem .surf files without extension')
    return parser.parse_args()

def get_initial_vertices(configfile):
    """Get config file options from section containing strings.
    :param config: ConfigParser object.
    :param section: Name of the section to be read.

    Args:
        configfile(str): The chosen config file.

    Returns:
        base1(int): The starting index for the mitral valve.
        base2(int): The starting index for the aortic valve.
        base3(int): The starting index for the pulmonary valve.
        base4(int): The starting index for the tricuspid valve.
        vertex_epi(int): The starting index for the epicardium.
        vertex_rv(int): The starting index for the right ventricle.
        vertex_lv(int): The starting index for the left ventricle.
    """
    config = configparser.ConfigParser()
    config.read(f"{configfile}")
    base1 = int(config['start_indices']['mv'])
    base2 = int(config['start_indices']['av'])
    base3 = int(config['start_indices']['pv'])
    base4 = int(config['start_indices']['tv'])
    vertex_epi = int(config['start_indices']['epi'])
    vertex_rv = int(config['start_indices']['rv'])
    vertex_lv = int(config['start_indices']['lv'])
    
    return base1, base2, base3, base4, vertex_epi, vertex_rv, vertex_lv



class EdgeToFaceMap:
    """Class containing methods to translate edges into faces.
    """
    def __init__(self):
        self.map = {}

    def construct_edge(self, node_a, node_b):
        """Inverts position of the nodes if the first one is bigger than the second one
        in order to construct an edge between them.
        
        Args: 
            node_a(int): The first node.
            node_b(int): The second node.

        Returns:
            node_a(int): The smaller node of the two.
            node_b(int): The bigger one.
        """
        if node_a > node_b:
            node_b, node_a = node_a, node_b
        return node_a, node_b

    def add_edge_if_missing(self, node_a, node_b, face_id):
        """Adds an edge between two nodes if it does not exist yet.
        
        Args:
            node_a(int): The first node.
            node_b(int): The second node.
            face_id(int): The id of the given face.
        
        Returns:
            None
        """
        map = self.map
        node_a, node_b = self.construct_edge(node_a, node_b)

        assert node_a < node_b
        edge = (node_a, node_b)
        if edge in map:
            map[edge].append(face_id)
        else:
            map[edge] = [face_id]

    def lookup(self, node_a, node_b):
        """Figures out if edge exists in the map.
        
        Args:
            node_a(int): The first node.
            node_b(int): The second node.
        """
        edge = self.construct_edge(node_a, node_b)
        if edge not in self.map:
            raise ValueError(f"Can't find edge ({node_a}, {node_b})")
        return self.map[edge]


class NormalCache(object):
    """Class for normals.
    """
    def __init__(self, faces, vertices):
        self.normals = {}
        self.faces = faces
        self.vertices = vertices

    def compute_normal(self, face_id):
        """Computes element normal.

        Args:
            face_id(int): The given face id.
        Returns:
            normal(np.ndarray): The element normal.
        """
        face = self.faces[face_id]
        l1 = self.vertices[face[0]] - self.vertices[face[1]]
        l2 = self.vertices[face[0]] - self.vertices[face[2]]
        
        normal = np.cross(l1, l2)
        normal /= np.linalg.norm(normal)

        return normal

    def lookup(self, face_id):
        """Looks up which face ids belong to the normal

        Args: 
            face_id(int): The given face id.
        """
        if face_id not in self.normals:
            normal = self.compute_normal(face_id)
            self.normals[face_id] = normal
            return normal
        return self.normals[face_id]
    
    
    

def grow_surface(source_face, face_graph : nx.graph, normal_cache: NormalCache, threshold):
    """Function that expands one vertex to the desired surface.
    
    Args:
        source_face(int64): The face id of a given vertex.
        face_graph(nx.classes.graph.Graph): A graph where the graph nodes are faces and graph edges are faces that share an edge.
        normal_cache(class): Class of normals.
        threshold(float): The given threshold.

    Returns:
        surface(set): Set containing all ids of the chosen surface.
    """

    surface = set()
    surface.add(source_face)
    faces_to_check = list()
    for neigh in face_graph.adj[source_face]:
        faces_to_check.append((source_face, neigh))

    while len(faces_to_check) > 0:
        # good_face is already on the surface, candidate_face should be checked
        good_face, candidate_face = faces_to_check.pop(0)
        if candidate_face in surface:
            continue

        good_normal = normal_cache.lookup(good_face)
        candidate_normal = normal_cache.lookup(candidate_face)
        
        dot_product = np.dot(good_normal, candidate_normal)
        if dot_product > 1:
            angle_rad = 0
        elif dot_product < -1:
            angle_rad = np.pi
        else:
            angle_rad = np.arccos(dot_product)
        
        if angle_rad < threshold:
            surface.add(candidate_face)
            # this candidate is also on the surface
            for neigh in face_graph.adj[candidate_face]:
                faces_to_check.append((candidate_face, neigh))

    return surface


if __name__ == "__main__":
    
    args = get_arguments(sys.argv)
    filename = args.input_filename
    
    # define name of config file
    # remove file extension and add config extension
    base_filename = filename[:-4] 
    config_file = filename[:-4]+".config"
    # get initial vertices from config file
    vertex_base1, vertex_base2, vertex_base3, vertex_base4, vertex_epi, vertex_rv, vertex_lv = get_initial_vertices(config_file)
    # Params
    base_depth_limit = 1 # identifies only surronding nodes around seed
    base_iterations = 20 # iterate multiple times to identify all connected nodes
    surface_threshold = 0.29 # 0.29 rad; 11 deg
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
        
    edge_to_face_map = EdgeToFaceMap()
    for i in range(0, num_faces):
        g.add_edge(int(faces[i][0]), int(faces[i][1]))
        g.add_edge(int(faces[i][1]), int(faces[i][2]))
        g.add_edge(int(faces[i][2]), int(faces[i][0]))
        g.nodes[int(faces[i][0])][i] = "key"
        g.nodes[int(faces[i][1])][i] = "key"
        g.nodes[int(faces[i][2])][i] = "key"

        edge_to_face_map.add_edge_if_missing(int(faces[i][0]), int(faces[i][1]), i)
        edge_to_face_map.add_edge_if_missing(int(faces[i][1]), int(faces[i][2]), i)
        edge_to_face_map.add_edge_if_missing(int(faces[i][2]), int(faces[i][0]), i)

    # construct a graph where the graph nodes are faces and graph edges are faces that share an edge
    face_graph = nx.Graph()
    for edge, face_ids in edge_to_face_map.map.items():
        assert len(face_ids) <= 2
        if len(face_ids) == 2:
            face_id_a, face_id_b = face_ids
            face_graph.add_edge(int(face_id_a), int(face_id_b))

    # Find face_id from vertex
    base1 = find_face_id(vertex_base1, faces, num_faces)
    base2 = find_face_id(vertex_base2, faces, num_faces)
    base3 = find_face_id(vertex_base3, faces, num_faces)
    base4 = find_face_id(vertex_base4, faces, num_faces)
    epi = find_face_id(vertex_epi, faces, num_faces)
    rv = find_face_id(vertex_rv, faces, num_faces)
    lv = find_face_id(vertex_lv, faces, num_faces)

    nx.bfs_tree(face_graph, source=base1[0], depth_limit=2).nodes()

    normal_cache = NormalCache(faces, vertices)
    normal_cache.lookup(base1[0])

    # Extract the surfaces for all bases
    base1_surface = grow_surface(base1[0], face_graph, normal_cache, surface_threshold)
    base2_surface = grow_surface(base2[0], face_graph, normal_cache, surface_threshold)
    base3_surface = grow_surface(base3[0], face_graph, normal_cache, surface_threshold)
    base4_surface = grow_surface(base4[0], face_graph, normal_cache, surface_threshold)

    base1_surface = list(base1_surface)
    base1_surface_array = np.array(base1_surface)
    base1_surface_array.sort()

    base2_surface = list(base2_surface)
    base2_surface_array = np.array(base2_surface)
    base2_surface_array.sort()

    base3_surface = list(base3_surface)
    base3_surface_array = np.array(base3_surface)
    base3_surface_array.sort()

    base4_surface = list(base4_surface)
    base4_surface_array = np.array(base4_surface)
    base4_surface_array.sort()
 
 
    face_tag = np.zeros(num_faces)
    face_tag_sum = face_tag
    tag = 1
    face_tag[base1_surface_array] = tag
    write_surfaces(f"{base_filename}_mv.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    face_tag_sum += face_tag

    face_tag = np.zeros(num_faces)
    tag = 1
    face_tag[base2_surface_array] = tag
    write_surfaces(f"{base_filename}_av.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    face_tag_sum += face_tag

    face_tag = np.zeros(num_faces)
    tag = 1
    face_tag[base3_surface_array] = tag
    write_surfaces(f"{base_filename}_pv.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    face_tag_sum += face_tag

    face_tag = np.zeros(num_faces)
    tag = 1
    face_tag[base4_surface_array] = tag
    write_surfaces(f"{base_filename}_tv.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    face_tag_sum += face_tag

    # Exclude bases from further processing
    face_tag = np.zeros(num_faces, dtype=int)
    for i in range(0, num_faces):
        if face_tag_sum[i] > 0:
            face_tag[i] = -1
    face_tag_exclbase = face_tag

    # Remove bases from graph
    # Tag exclusion is not enough since the graph is still connected
    # Graph nodes must be removed to get a segmented graph for further BFS iterations
    for i in range(0, num_faces):
        if face_tag[i] == -1:
            try:
                g.remove_node(int(faces[i][0]))
                g.remove_node(int(faces[i][1]))
                g.remove_node(int(faces[i][2]))
            except Exception:
                pass
    print("[info] extracting epicardium")
    face_tag[epi[0]] = 1
    face_tag, tag = apply_bfs_algorithm(face_tag, faces, num_faces, g, surface_depth_limit)
    write_surfaces(f"{base_filename}_epi.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)

    # Exclude epi from further processing
    for i in range(0, num_faces):
        if face_tag[i] > 0:
            face_tag[i] = -1

    print("[info] extracting left ventricular endocardium")
    face_tag[lv[0]] = 1
    face_tag, tag = apply_bfs_algorithm(face_tag, faces, num_faces, g, surface_depth_limit)
    write_surfaces(f"{base_filename}_endo_lv.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)
    
    # Exclude lv from further processing
    for i in range(0, num_faces):
        if face_tag[i] > 0:
            face_tag[i] = -1

    print("[info] extracting right ventricular endocardium")
    face_tag[rv[0]] = 1
    face_tag, tag = apply_bfs_algorithm(face_tag, faces, num_faces, g, surface_depth_limit)
    write_surfaces(f"{base_filename}_endo_rv.ply", tag, face_tag, faces, num_faces, vertices, num_vertices)

   