import numpy as np

def louvain_community_detection(adj_list):    
    def adj_list_to_matrix(adj_list):
    
        size = max(max(nodes) for nodes in adj_list) + 1
        adj_matrix = np.zeros((size, size), dtype=int)
        for node, neighbors in enumerate(adj_list):
            for neighbor in neighbors:
                adj_matrix[node][neighbor] = 1
        return adj_matrix
    def modularity_gain(node, community, adj_matrix, degrees, m):
        community_nodes = np.where(community == community[node])[0]
        in_degree_sum = adj_matrix[node, community_nodes].sum()
        degree_sum = degrees[community_nodes].sum()
        node_degree = degrees[node]
        
        gain = (in_degree_sum /(2* m)) - (node_degree * degree_sum) / (2 * m**2)
        return gain

    def phase_one(communities, adj_matrix, degrees, m):
        improved = True
        while improved:
            improved = False
            for node in range(n):
                best_community = communities[node]
                best_gain = 0
                original_community = communities[node]
                
                communities[node] = -1
                
                for neighbor in np.where(adj_matrix[node] > 0)[0]:
                    communities[node] = communities[neighbor]
                    gain = modularity_gain(node, communities, adj_matrix, degrees, m)
                    
                    if gain > best_gain:
                        best_gain = gain
                        best_community = communities[neighbor]
                
                communities[node] = best_community
                
                if best_community != original_community:
                    improved = True
        print(f"communities after one iteration of louvain algo are {len(np.unique(communities))}")        
        return communities

    def phase_two(communities, adj_matrix):
        """Second phase of Louvain, aggregating nodes into their communities."""
        unique_communities = np.unique(communities)
        new_adj_matrix = np.zeros((len(unique_communities), len(unique_communities)))
        new_node_map = {old: new for new, old in enumerate(unique_communities)}
        
        for i in range(len(adj_matrix)):
            for j in range(len(adj_matrix)):
                if communities[i] == communities[j]:
                    new_i = new_node_map[communities[i]]
                    new_j = new_node_map[communities[j]]
                    new_adj_matrix[new_i, new_j] += adj_matrix[i, j]
        
        return new_adj_matrix
    adj_matrix=adj_list_to_matrix(adj_list)
    n = len(adj_matrix)
    adj_matrix=np.array(adj_matrix)
    
    degrees = np.sum(adj_matrix,axis=1)
    m = degrees.sum()/2
    communities = np.arange(n)

        
    communities = phase_one(communities, adj_matrix, degrees, m)
    new_adj_matrix = phase_two(communities, adj_matrix)
    return communities,new_adj_matrix
    # Initial setup
    


