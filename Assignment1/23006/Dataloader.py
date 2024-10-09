import numpy as np
from collections import defaultdict


def import_wiki_vote_data(file_path):
    edges=[]
    with open(file_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            i_,j_=map(int, line.split())
            i=max(i_,j_)
            j=min(i_,j_)
            edges.append([i,j])
        edges_array =np.array(edges)
        
        edges_array=np.unique(np.sort(edges_array, axis=1), axis=0 )
        
        return edges_array
    
def mapping(edges_array):
    
    nodes=np.unique(edges_array)
    node_mapping = {}
    for idx, node in enumerate(nodes):
        node_mapping[node] = idx
    remapped_edges = [(node_mapping[i], node_mapping[j]) for i, j in edges_array]
    return remapped_edges, node_mapping

def create_reverse_mapping(node_mapping):
    reverse_mapping = {v: k for k, v in node_mapping.items()}
    return reverse_mapping
    
def edge_list_to_adjacency_list(remapped_edges):
    adjacency_list = [[] for _ in range(len(np.unique(remapped_edges)))]

    for edge in remapped_edges:
        i,j=edge
        adjacency_list[i].append(j)
        adjacency_list[j].append(i)
    
    return adjacency_list

edge_list = import_wiki_vote_data("../data/Wiki-Vote.txt")[:50]
remapped, mapped = mapping(edge_list)
final_list = edge_list_to_adjacency_list(remapped)



def import_lastfm_asia_data(file_path):
    # Load the CSV file into a numpy array
    data = np.array(np.genfromtxt(file_path, delimiter=',', dtype=int)[1:])
    return data
