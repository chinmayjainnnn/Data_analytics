from collections import defaultdict, deque
from Dataloader import final_list
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from time import time 
from scipy.cluster.hierarchy import linkage, dendrogram

def bfs_shortest_paths(graph, source):
    """ 
    Perform BFS to find the shortest paths from the source to all other nodes.
    Returns:
        sigma: Number of shortest paths from source to each node.
        delta: Dependency of nodes on the source.
        predecessors: Predecessor nodes in the shortest paths.
    """
    sigma = defaultdict(int)
    delta = defaultdict(float)
    predecessors = defaultdict(list)
    sigma[source] = 1
    
    queue = deque([source])
    distances = {source: 0}
    visited = set([source])
    
    while queue:
        v = queue.popleft()
        for w in graph[v]:
            if w not in visited:
                queue.append(w)
                visited.add(w)
                distances[w] = distances[v] + 1
            if distances[w] == distances[v] + 1:
                sigma[w] += sigma[v]
                predecessors[w].append(v)
    
    return sigma, delta, predecessors, distances

def calculate_betweenness_mp(graph, start, end):
    betweenness = defaultdict(float)
    
    for node in range(start, end):
        sigma, delta, predecessors, distances = bfs_shortest_paths(graph, node)
        stack = sorted(distances, key=distances.get, reverse=True)
        
        for w in stack:
            for v in predecessors[w]:
                c = (sigma[v] / sigma[w]) * (1 + delta[w])
                betweenness[(v, w)] += c
                betweenness[(w, v)] += c
                delta[v] += c
    # for edge in betweenness:
    #     print(f"edge={edge}")
    #     break
    for edge in betweenness:
        betweenness[edge] /= 2.0
    
    combined_betweenness = defaultdict(float)
    for (u, v), between in betweenness.items():
        # Sort the nodes to ensure that (u, v) and (v, u) are treated the same
        sorted_edge = tuple(sorted((u, v)))
        combined_betweenness[sorted_edge] += between


    return combined_betweenness


def calculate_betweenness(graph):
    """
    Calculate the betweenness centrality for all edges.
    Returns:
        betweenness: A dictionary with edges as keys and betweenness centrality as values.
    """
    betweenness_final_list = []
    start_time = time()
    betweenness = defaultdict(float)
    n = len(graph)
    
    num_processes = 80
    # Divide the range of nodes among processes
    chunk_size = (n + num_processes - 1) // num_processes  # To handle uneven division
    ranges = [(i, min(i + chunk_size, n)) for i in range(0, n, chunk_size)]

    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(calculate_betweenness_mp, [(graph, start, end) for start, end in ranges])

    for local_betweenness in results:
        for edge, value in local_betweenness.items():
            betweenness[edge] += value

    print(f"Time: {time()-start_time}")
    return betweenness

def bfs(start, visited, adjacency_list):
    queue = deque([start])
    connected_component = []
    
    while queue:
        node = queue.popleft()
        if not visited[node]:
            visited[node] = True
            connected_component.append(node)
            for neighbor in adjacency_list[node]:
                if not visited[neighbor]:
                    queue.append(neighbor)
    
    return connected_component

def find_components(adjacency_list):
    visited = [False] * len(adjacency_list)
    connected_components = []
    
    for node in range(len(adjacency_list)):
        if not visited[node]:
            connected_component = bfs(node, visited, adjacency_list)
            connected_components.append(connected_component)
    
    return connected_components

def girvan_newman(adj_list):
    graph = defaultdict(set)
    for i, neighbors in enumerate(adj_list):
        for neighbor in neighbors:
            graph[i].add(neighbor)
            graph[neighbor].add(i)
    
    partition = []
    while any(graph.values()):
        community_map,graph = girvan_newman_one_level(graph)
        partition.append(community_map)
        np.save('part1.npy',np.array(partition))
    partition = np.concatenate(partition,axis=1)
    return partition


def girvan_newman_one_level(graph):
    components = find_components(graph)
    old_comp = len(components)
    while True:
        betweenness = calculate_betweenness(graph)
        # print(f"{betweenness=}")
        if not betweenness:
            break
        
        # Remove the edge with the highest betweenness centrality
        max_edge = max(betweenness, key=betweenness.get)
        graph[max_edge[0]].remove(max_edge[1])
        graph[max_edge[1]].remove(max_edge[0])
        print(max_edge[0],max_edge[1])
        components = find_components(graph)
        print(len(components))
        if len(components) == old_comp + 1:
            break

    community_map = np.zeros(len(graph), dtype=int)
    for component in components:
        min_node = min(component)
        for node in component:
            community_map[node] = min_node
    
    return community_map.reshape(-1,1),graph


def plot_dendrogram(Z, filename='dendrogram1.png'):
    
    plt.figure(figsize=(10, 7))
    Z = linkage(Z, 'ward')
    dendrogram(Z)
    plt.title('Dendrogram for Communities')
    plt.xlabel('Node Index')
    plt.ylabel('Distance')
    
    # Save the dendrogram as a .png file
    plt.savefig(filename)
    plt.show()