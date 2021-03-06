import networkx as nx
import argparse
import sys
import pandas as pd
import os
import numpy as np
import csv
import itertools
from functools import reduce
import random

# Read network

def read_network(network):
    
    print('\n')
    print ("** Constructing network from  file...")
    with open(network, 'r') as network_file:
        G=nx.read_edgelist(network_file, delimiter='\t', nodetype=int)
        G.remove_edges_from(nx.selfloop_edges(G))
    print ("\tNetwork contains %s nodes and %s edges\n" %(G.number_of_nodes(), G.number_of_edges()))
    return G


# Read seed nodes file (seed nodes are the subset of nodes in the network [i.e, identified louvain community]

def read_seed(ifile):
    
    print("** Parsing seed nodes...")
    seed_nodes = set()
    with open(ifile, 'r') as seed_file:
        reader = csv.reader(seed_file)
        for seed in reader:
            seed1 = int(seed[0])
            seed_nodes.add(seed1)
    print ("\tNo of seed nodes: %s\n" %(len(seed_nodes)))
    return seed_nodes


# Generate subgraph with seed nodes

def gen_subgraph(G, seeds):
    
    g = nx.Graph(G.subgraph(seeds))
    return g


# Calculated degree centrality of the network

def deg_centrality(g):

    dc = nx.degree_centrality(g)
    DC = pd.DataFrame(dc.items(), columns=['node', 'DC'], dtype=int)
    #print (DC)
    return DC

def avg_deg_centrality(DC):
    
    DC_mean = DC['DC'].mean()
    #print (DC_mean)
    return DC_mean


# Calculated closeness centrality of the network

def close_centrality(g):

    cc = nx.closeness_centrality(g)
    CC = pd.DataFrame(cc.items(), columns=['node', 'CC'], dtype=int)
    #print (CC)
    return CC

def avg_close_centrality(CC):

    CC_mean = CC['CC'].mean()
    #print (CC_mean)
    return CC_mean

# Calculated eigenvector centrality of the network

def eigen_centrality(g):

    ec = nx.eigenvector_centrality_numpy(g)
    EC = pd.DataFrame(ec.items(), columns=['node', 'EC'], dtype=int)
    #print (EC)
    return EC

def avg_eigen_centrality(EC):

    EC_mean = EC['EC'].mean()
    #print (EC_mean)
    return EC_mean


# Calculate average centralities for random networks

def random_comparisons(G, seeds, n, DC_mean, CC_mean, EC_mean):

    print("** Performing random simulations to test the significance of observed centrality scores.")
    all_genes = G.nodes()
    len_of_seed = len(seeds)
    
    dc_rand = []
    cc_rand = []
    ec_rand = []

    for i in range(1,n+1):
     
        rand = set(random.sample(all_genes,len_of_seed))
        rg = gen_subgraph(G, rand)

        rand_dc = deg_centrality(rg)
        #rand_dc_df = pd.DataFrame(rand_dc.items(), columns=['node', 'DC'], dtype=int)
        avg_rand_dc = avg_deg_centrality(rand_dc)
        dc_rand.append(avg_rand_dc)

        rand_cc = close_centrality(rg)
        #rand_cc_df = pd.DataFrame(rand_cc.items(), columns=['node', 'CC'], dtype=int)
        avg_rand_cc = avg_close_centrality(rand_cc)
        cc_rand.append(avg_rand_cc)

        rand_ec = eigen_centrality(rg)
        #rand_ec_df = pd.DataFrame(rand_ec.items(), columns=['node', 'EC'], dtype=int)
        avg_rand_ec = avg_eigen_centrality(rand_ec)
        ec_rand.append(avg_rand_ec)
    
    # calculate p-values

    dc_true = []
    cc_true = []
    ec_true = []

    for dc in dc_rand:
        if dc >= DC_mean:
            dc_true.append(1)
        else:
            pass
    
    for cc in cc_rand:
        if cc >= CC_mean:
            cc_true.append(1)
        else:
            pass

    for ec in ec_rand:
        if ec >= EC_mean:
            ec_true.append(1)
        else:
            pass
    
    dct_count = dc_true.count(1)
    cct_count = cc_true.count(1)
    ect_count = ec_true.count(1)

    pval_dc = dct_count/n
    pval_cc = cct_count/n
    pval_ec = ect_count/n

    pvalues = [pval_dc, pval_cc, pval_ec]

    print("\tDone.")
    
    return pvalues


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-n", "--network", required=True, help="Network file. A tab-separated file containing a list of numeric ids of two interacting proteins (REQUIRED)", dest = "network")
    parser.add_argument("-g", "--genes", required=True, help="The significant_eqtls.txt file or a txt file containing a list of gene entrez ids (REQUIRED)", dest= "seeds")
    parser.add_argument("-r", "--random", required=False, default = 1000, help="Number of random simulations to be calculated (default 1000)", dest= "randm")
    parser.add_argument("-o", "--output", required=True, help="directory to write results (REQUIRED)", dest = "out")
    global args
    args = parser.parse_args()

    if args.seeds == 'none':
        message = ''' 
        ERROR: you must specify input file with gene ids.
        For more details: use measure_centrality.py -h
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
        For more details: use find_significant_module.py -h
        '''
        print (message)
        sys.exit(0)

    # Load network file

    G = read_network(args.network)

    # Read seed nodes file

    seeds = read_seed(args.seeds)

    # Generate subgraph with seed genes 
    
    print("** Generating subgraph from seed nodes...")
    g = gen_subgraph(G, seeds)
    print("\tSubgraph generated.")
    print("\tSubgraph contains %s nodes and %s edges." %(g.number_of_nodes(), g.number_of_edges()))
    print("\tWriting edgelist of subgraph to output file (cluster_subgraph.txt).\n")
    
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    
    path = os.path.join(args.out, 'cluster_subgraph.txt')
    f = open(path, 'wb')
    nx.write_edgelist(g, path = f, delimiter='\t', data=False, encoding='utf-8')

    # Calculate degree centrality

    DC = deg_centrality(g)

    DC_mean = avg_deg_centrality(DC)

    # Calculate closeness centrality

    CC = close_centrality(g)

    CC_mean = avg_close_centrality(CC)

    # Calculate eigenvector centrality

    EC = eigen_centrality(g)

    EC_mean = avg_eigen_centrality(EC)

    # Combine all centrality values into a single dataframe
    
    C = [DC, CC, EC]

    scores = reduce(lambda left,right: pd.merge(left,right, on = ['node'], how = 'inner'), C)
    scores = scores.apply(lambda x: round(x, 3))
    outfile = os.path.join(args.out, 'centrality_scores.txt')
    scores.to_csv(outfile, sep = "\t", index=False, header=True)

    # calculate centralities on random seeds

    pval = random_comparisons(G, seeds, int(args.randm), DC_mean, CC_mean, EC_mean)
    
    pval_outfile = os.path.join(args.out, 'pvalues.txt')
    pval_open = open(pval_outfile, 'w')
    writer = csv.writer(pval_open, delimiter = '\t')
    header = ['DC', 'CC', 'EC']
    writer.writerow(header)
    writer.writerows([pval])

