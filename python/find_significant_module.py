import networkx as nx
import argparse
import sys
import csv 
import pandas as pd
import os
import community as community_louvain
import qstest as qs

#######################################################################################################
'''
    * Read in network file as edgelist
    * Read network from an external file
    * Network file format: Must contain only two columns with a list of Entrez gene IDs representing 
      interacting proteins
'''
#######################################################################################################


def read_network(network):
    print('\n')
    print ("** Loading and processing network file...")
    with open(network, 'r') as network_file:
        G=nx.read_edgelist(network_file, delimiter='\t', nodetype=int)
    print ("\tNetwork file processing done")
    print ("\tNetwork contains %s nodes and %s edges\n" %(G.number_of_nodes(), G.number_of_edges()))
    return G

#######################################################################################################
'''
    * To remove self-links
'''    
#######################################################################################################


def remove_self_links(G):
    G.remove_edges_from(nx.selfloop_edges(G))


#######################################################################################################
'''
    * Input file should be a tab-delimited file
    * First column of the input gene file should contain Entrez IDs or any integer as gene identifier
    * Input gene list can have multiple columns but only first column will be considered
'''
#######################################################################################################


def input_gene_list(file_in):
    print('** Reading input gene list...')
    gene_list = set()
    with open(file_in, 'r') as gene_file:
        reader = csv.reader(gene_file, delimiter = "\t")
        for row in reader:
            row1 = int(row[0])
            #gene_df.join(row[0].unique())
            gene_list.add(row1)
    print ("\tNo of genes in the input gene list: %s" %(len(gene_list)))
    return gene_list


#######################################################################################################
'''
    * Filter out genes from the input file that are not present in the network
'''
#######################################################################################################


def filter_input_genes(gene_list, G):
    all_genes_in_network = set(G.nodes())
    input_genes_in_network = set(gene_list & all_genes_in_network)
    return input_genes_in_network


#######################################################################################################
'''
    * Apply louvain algorithm to find communities within the network contructed
    * Communities identified by the algorithm with only one node will be removed from further analysis
'''
#######################################################################################################


def louvain(g):
    
    coms = community_louvain.best_partition(g) 
    # coms return dictionaries of nodes and communities as key and value, eg:{'MTCH2': 0, 'PPIL1': 1,}
    comm = pd.DataFrame(coms.items(), columns=['node', 'cluster'], dtype=int)
    groups = comm.groupby('cluster')['node'].apply(list)
    
    communities = []
    for grp in groups:
        if len(grp)>1:
            communities.append(grp)
    return communities


#######################################################################################################
'''
    * Checks the quality of communities
    * Calculate scores for all the communities identified in the above step
    * Used by qstest as quality function to check statistical significance of the communitites
'''
#######################################################################################################


def find_modularity(g, communities):

    deg = g.degree(communities)
    q = 0
    D = 0
    for i in communities:
        for j in communities:
            if g.has_edge(i, j) == False:
                continue
            q += 1.0
        D += deg[i]

    M = g.size() / 2 # total no of edges /2
    q = (q - D * D / (2.0 * M)) / (2 * M)
    return q

#######################################################################################################    
'''
    * Calculate size of the communities
    * Used by qstest as size function to check statistical significance of the communitites
'''
#######################################################################################################


def find_size(g, communities):
    s = sum([x[1] for x in g.degree(communities)])
    return s


#######################################################################################################
''' 
    * Write all communitites to the file "tested_communities.txt"
    * Output file will have four columns:
          column 1 : list of gene ids in a community
          column 2 : number of genes in a community
          column 3 : True/False denotes the significance
          column 4 : p-value of the communities
'''
#######################################################################################################


def write_sig_comm_tofile(sg, pval, communities, out):

    path = os.path.join(out, 'tested_communities.txt')
    output = open(path, 'w')
    
    size = []
    for c in communities:
        size.append(len(c))

    for i,j,k,l in zip(communities,size,sg,pval):
        output.write(str(i)+"\t"+str(j)+"\t"+str(k)+"\t"+str(l)+"\n")
    output.close()


#######################################################################################################
''' * End of functions '''
#######################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--genes", required=True, help="A text file containing a list of gene entrez ids (REQUIRED)", dest= "genes")
    parser.add_argument("-n", "--network", required=True, help="Network file. A tab-separated file containing a list of numeric ids of two interacting proteins (REQUIRED)", dest = "network")
    parser.add_argument("-o", "--output", required=True, help="directory to write results (REQUIRED)", dest = "out")
    parser.add_argument("-r", "--random_simulations", required=False, default = 1000, help="No of random simulations to run (OPTIONAL)", dest = "random")
    parser.add_argument("-t", "--no_of_threads", required=False, default = 12, help="No of threads to use (OPTIONAL)", dest = "threads")

    global args
    args = parser.parse_args()

    
    if args.genes == 'none':
        message = ''' 
        ERROR: you must specify input file with gene set.
        For more details: use find_significant_module.py -h 
        '''
        print (message)
        sys.exit(0)

    
    if args.network == 'none':
        message = ''' 
        ERROR: you must specify a network file.
        For more details: use find_significant_module.py -h 
        '''
        print (message)
        sys.exit(0)
    
    
    if args.out == 'none':
        message = ''' 
        ERROR: you must specify path to save the output file.
        For more details: use construct.py -h 
        '''
        print (message)
        sys.exit(0)

    
    # Load network file
    G = read_network(args.network)
    remove_self_links(G) 
        
    # Load input gene file
    all_input_genes = input_gene_list(args.genes) 

    
    # Remove genes in the input file that are not in the network
    input_genes_in_network = filter_input_genes(all_input_genes, G) 
    print ("\tUsing only %s genes from input gene list that are found in the network" %(len(input_genes_in_network)))
    
    
    # Generating subgraph with genes from the input file and write it to the output file as edgelist
    
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    
    g = nx.Graph(G.subgraph(input_genes_in_network))
    if g.number_of_edges() != 0:
        print ("\tGenerating subgraph from the input genes and writing to subgraph.txt in the output directory\n") 
        path = os.path.join(args.out, 'subgraph.txt')
        f = open(path, 'wb')
        nx.write_edgelist(g, path = f, delimiter='\t', data=False, encoding='utf-8')

    else:
        print ("        Subgraph has no edges.")

    
    # Identify communities
    print ("** Finding communities (disease modules)...")
    communities = louvain(g)
    coms = community_louvain.best_partition(g)
    C = max(coms.values()) + 1
    print ("        Done.")

    # Check statistical significance of the communities over randomly generated communities
    print ("\n** Testing significance of the communities....")


    stat_test = qs.qstest(g, communities, find_modularity, find_size, louvain, num_of_rand_net = int(args.random), num_of_thread = int(args.threads))
    values = [list(x) for x in stat_test]
    sg = values[0]
    pval = values[1]


    # Write significant communities to file
    out = write_sig_comm_tofile(sg, pval, communities, args.out)
    

    # Create logfile with miscellaneous details
    ofile = os.path.join(args.out, 'find_cluster_logfile.txt')
    outfile = open(ofile, 'w')
    outfile.write(" ** Network Details:\n")
    outfile.write("         Network contains %s nodes and %s edges\n" %(G.number_of_nodes(), G.number_of_edges()))
    outfile.write(" ** Input gene file:\n")
    outfile.write("         No of genes in the input gene list: %s\n" %(len(all_input_genes)))
    outfile.write("         Using only %s genes from input gene list that are found in the network\n" %(len(input_genes_in_network)))
    outfile.write("         %s clusters were identified\n" %(C))
    outfile.write("         subgraph generated with input genes is saved to subgraph.txt file if the subgraph has atleast one edge\n")
    outfile.write(" ** Random simulation results:\n")
    outfile.write("         Significance of the clusters (boolean): {}\n".format(sg)) 
    outfile.write("         pvalues corresponding to the clusters: {}\n".format(pval))
    outfile.close()

    print ("        Done.\n")
