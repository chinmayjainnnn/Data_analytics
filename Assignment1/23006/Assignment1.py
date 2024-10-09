import numpy as np
from Dataloader import *
from Girvan_Newman import *
import sys
from louvain import *

sys.setrecursionlimit(1000000000)




if __name__ == "__main__":
    

    nodes_connectivity_list_wiki = import_wiki_vote_data("../data/Wiki-Vote.txt")
    remapped, mapped = mapping(nodes_connectivity_list_wiki)
    final_list = edge_list_to_adjacency_list(remapped)
    community_mat_wiki = girvan_newman(final_list)
    plot_dendrogram(community_mat_wiki)
        

    graph_partition_wiki  = girvan_newman_one_level(final_list)

    
    nodes_connectivity_list_lastfm = import_lastfm_asia_data("../data/lastfm_asia_edges.csv")
    final_list = edge_list_to_adjacency_list(nodes_connectivity_list_lastfm)
    graph_partition_louvain_lastfm = louvain_community_detection(final_list)
    graph_partition_lastfm = girvan_newman_one_level(final_list)
    community_mat_lastfm = girvan_newman(nodes_connectivity_list_lastfm)
    plot_dendrogram(community_mat_lastfm)


    nodes_connectivity_list_wiki = import_wiki_vote_data("../data/Wiki-Vote.txt")
    remapped, mapped = mapping(nodes_connectivity_list_wiki)
    final_list = edge_list_to_adjacency_list(remapped)
    z=np.load('/raid/home/chinmayjain/NLP/Data_analytics/code/part.npy')

    z=np.concatenate(z,axis=1)
    stopping_criteria(z,final_list)
####################################################################################################
    
    # z=np.load('/raid/home/chinmayjain/NLP/Data_analytics/code/part.npy')
    # z=np.concatenate(z,axis=1)
    # plot_dendrogram(z)





########################################################################################################




    ############ Answer qn 1-4 for wiki-vote data #################################################
    # Import wiki-vote.txt
    # nodes_connectivity_list is a nx2 numpy array, where every row 
    # is an edge connecting i->j (entry in the first column is node i, 
    # entry in the second column is node j)
    # Each row represents a unique edge. Hence, any repetitions in data must be cleaned away.
    

    # This is for question no. 1
    # graph_partition: graph_partitition is a nx1 numpy array where the rows corresponds to nodes in the network 
    #                   (0 to n-1) and  the elements of the array are the community ids of the corressponding nodes.
    #                  Follow the convention that the community id is equal to the lowest nodeID in that community.
    

    # This is for question no. 2. Use the function 
    # written for question no.1 iteratetively within this function.
    # community_mat is a n x m matrix, where m is the number of levels of Girvan-Newmann algorithm and n is the number of nodes in the network.
    # Columns of the matrix corresponds to the graph_partition which is a nx1 numpy array, as before, corresponding to each level of the algorithm. 
    
    #///////////////////////////////
    
    # This is for question no. 3
    # Visualise dendogram for the communities obtained in question no. 2.
    # Save the dendogram as a .png file in the current directory.
    

    # This is for question no. 4
    # run one iteration of louvain algorithm and return the resulting graph_partition. The description of
    # graph_partition vector is as before. Show the resulting communities after one iteration of the algorithm.
    # graph_partition_louvain_wiki = louvain_one_iter(nodes_connectivity_list_wiki)


    ############ Answer qn 1-4 for bitcoin data #################################################
    # Import lastfm_asia_edges.csv
    # nodes_connectivity_list is a nx2 numpy array, where every row 
    # is an edge connecting i<->j (entry in the first column is node i, 
    # entry in the second column is node j)
    # Each row represents a unique edge. Hence, any repetitions in data must be cleaned away.