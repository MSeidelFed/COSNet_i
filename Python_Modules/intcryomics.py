#!/usr/env/bin/python3
"""
[sample_network.py]: 
Usage: python3 intcryomics.py [filewithedgelist] [sigfile] [walklength] [iterationnum]
"""
## IMPORTS
from sys import argv
from scipy import stats
import operator
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import collections
import os
import random
import glob
import math
import pickle
import itertools
from statsmodels.sandbox.stats.multicomp import multipletests
import argparse
## FUNCTIONS
def get_network_from_edges(edgelistfile):
    """ """
    G = nx.read_edgelist(edgelistfile, data=(('weight',float),))
    print(f'List of nodes: {G.nodes()}\n')
    #for node in G.nodes():
        #print(f'Neighbors of node: {node}')
        #print(list(nx.all_neighbors(G,node)))
    adjmtx = nx.adj_matrix(G)
    adjmtx = adjmtx.todense()
    adjmtx = np.array(adjmtx, dtype = int)
    #degmtx = np.diag(np.sum(adjmtx, axis=0))
    #np.savetxt('degmatrix.dat', degmtx, delimiter=',')
    np.savetxt('adjmatrix.dat', adjmtx, delimiter=',') 

    return G, adjmtx

def transit_design(adjacent_mtx):
    """
    Returns matrix of transit probabilities for the random walk
    where probabilities are calculated based on the network weighting scheme:
        P = w_{x,y}/w_{x} where w_{x} is the sum of all weights for outgoing edges from node x
    """
    summtx = np.diag(np.sum(adjacent_mtx, axis=0))
    transmtx = np.dot(np.linalg.inv(summtx),adjacent_mtx)
    np.savetxt('Tmatrix.dat', transmtx, delimiter = ',')    
    return transmtx

def random_walk(transit_mtx, startnodelist, nodelist, walklength):
    """ 
    Implements random walk based on random starts (if given a multi-element list)
    or a particular startnode (list with one element), for a certain walklength,
    based on probabilities as given in the transit matrix
    """
    startnode = random.choice(startnodelist)
    startnode_idx = nodelist.index(startnode)
    v_i = np.zeros(len(nodelist))
    v_i[startnode_idx] = 1
    visited = [startnode_idx]
    visitednames = [startnode]
    for i in range(walklength):
        v_i = np.dot(v_i, transit_mtx)
        m = max(v_i)
        nextnode_idx = random.choice([i for i, j in enumerate(v_i) if j == m])
        visited.append(nextnode_idx)
        visitednames.append(nodelist[nextnode_idx])
        v_i = np.zeros(len(nodelist))
        v_i[nextnode_idx] = 1
    return visited, visitednames

def generateSequence(startIndex, transitionMatrix, path_length, alpha = 0.1):
    """ """
    result = [startIndex]
    current = startIndex

    for i in range(0, path_length):
        if random.random() < alpha:
            nextIndex = startIndex
        else:
            probs = transitionMatrix[current]
            nextIndex = np.random.choice(len(probs), 1, p=probs)[0]

        result.append(nextIndex)
        current = nextIndex

    return result

def random_walk_new(G, transitmtx, iter_num, walk_length):
    """ """
    nodes = list(G.nodes())
    WalksDict = {}
    count = 0

    for i in range(0, len(nodes)):
        for j in range(0, iter_num):
            indexList = generateSequence(i, transitmtx, walk_length)
            entryname = (nodes[i], count)
            WalksDict[entryname] = [nodes[tmp] for tmp in indexList]
            count += 1

    return WalksDict

def get_significance_values(filewithsigs):
    """ 
    Reads in a two-column file with proteins and binary encoding
    for significant or non-significant abundance changes
    
    Returns a dictionary of form: { prot : 0/1 }
    """
    sig_dict = {}
    with open(filewithsigs, 'r') as f:
        for line in f:
            line_elems = line.strip().split(' ')
            sig_dict[line_elems[0]] = int(line_elems[1])
    return sig_dict

def get_sig_vals_full(filewithsigs):
    """
    Gives back form of {uL1:[0,0,0],uL2:[1,1,1]...}
    """
    sig_dict_full = {}
    with open(filewithsigs,'r') as f:
        for line in f:
            line_elems = line.strip().split(' ')
            if not line_elems[0] in sig_dict_full:
                sig_dict_full[line_elems[0]] = [int(line_elems[1])]
            else:
                sig_dict_full[line_elems[0]].append(int(line_elems[1]))
    return sig_dict_full
                

def iterate_random_walk(nodelist, iteratenumber, transitmtx, walklength):
    """
    Iteratively runs the random walk based on a given iteration number
    """
    WalksDict = {}
    count = 0
    for j in range(len(nodelist)):
        startnode = [nodelist[j]]
        for i in range(iteratenumber):
            visit_list, visitednames = random_walk(transitmtx, startnode, nodelist, walklength)
            entryname = (startnode[0], count)
            WalksDict[entryname] = visit_list
            print(f'{entryname} {visit_list}')
            count += 1
    return WalksDict

def binom_test_regions(dictofregions, Nodelist, SigDict, siglvl=0.05):
    """ 
    Tests each region based on the binomial test with p = 18/80 = 0.225
    """
    TestResDict = {}
    for ID, region in dictofregions.items():
        sigcount = 0
        for protid in region:
            prot = Nodelist[protid]
            if int(SigDict[prot]) == 1:
                sigcount += 1
        [M, n, N] = [len(region), sigcount, math.ceil(len(region)/2)]
        hpd = stats.hypergeom(M,n,N)
        x = np.arange(0, n+1)
        pmf_prots = hpd.pmf(x)
        print(pmf_prots)
        print(stats.hypergeom.sf(math.ceil(len(region)/4), M, n, N))
        pval_binom = stats.binom_test(sigcount, n=len(region), p=0.225)
        TestResDict[ID] = pval_binom
    return TestResDict

def fishers_exact_test_regions(dictofregions, Nodelist, SigDict):
    """
    have to count how many actually in one region (including rps)
    and also how many sig
    """
    TestResDict = {}
    PvalList = []
    #totnumsig = sum(SigDict.values())
    masterlist = list(itertools.chain.from_iterable(SigDict.values()))
    totnumsig = sum(masterlist)
    print(f'Total number significant: {totnumsig}')
    for ID, region in dictofregions.items():
        sigcount = 0
        regcount = 0
        for protid in region:
            prot = Nodelist[protid]
            regcount += len(SigDict[prot])
            for sigval in SigDict[prot]:
                if sigval == 1:
                    sigcount += 1
            #if int(SigDict[prot]) == 1:
            #    sigcount += 1
        #not_sig_in_reg = len(region) - sigcount
        #sig_not_in_reg = totnumsig - sigcount
        #not_sig_not_in_reg = len(Nodelist) - len(region) - sig_not_in_reg
        #oddsratio, pval = stats.fisher_exact([[sigcount,sig_not_in_reg],[not_sig_in_reg, not_sig_not_in_reg]])

        not_sig_in_reg = regcount - sigcount
        sig_not_in_reg = totnumsig - sigcount
        not_sig_not_in_reg = len(masterlist) - regcount - sig_not_in_reg
        oddsratio, pval = stats.fisher_exact([[sigcount,sig_not_in_reg],[not_sig_in_reg, not_sig_not_in_reg]])
        TestResDict[ID] = (oddsratio, pval)
        PvalList.append(pval)

    return TestResDict, PvalList

def get_minimum_set_cover(nodelist, listofsubsets):
    """
    Implements the minimum set cover algorithm to find non-overlapping sets
    out of the 80 ribosomal sampled regions
    """
    indices = set(range(len(nodelist)))
    elems = set(e for s in listofsubsets for e in s)
    if elems != indices:
        return None
    covered = set()
    cover = []
    while covered != elems:
        subset = max(listofsubsets, key=lambda s: len(s - covered))
        cover.append(subset)
        covered |= subset
    return cover

def remove_overlaps(coverlist, dictofoverlaps):
    """
    Reassigns overlapping regions to smaller regions
    """
    dictofcover = {}
    for i in range(len(coverlist)):
        dictofcover[i] = list(coverlist[i])

    update_ov_dict = dictofoverlaps
    #print('Original coverlist')
    #print(coverlist)
    #print('Original overlaps')
    #print(dictofoverlaps)

    for i in range(len(dictofoverlaps)):
        if not update_ov_dict: #check if there are any overlaps left
            break
        else:
            pairs = list(update_ov_dict.keys())
            pair_to_fix = random.choice(pairs)
            overlap = update_ov_dict[pair_to_fix]
            #print(f'Overlap between {pair_to_fix}: {overlap}')
            regio1 = dictofcover[pair_to_fix[0]]
            regio2 = dictofcover[pair_to_fix[1]]
 
            if len(regio1) > len(regio2):
                for prot in overlap:
                    if prot in regio1:
                        regio1.remove(prot)
                    if not prot in regio2:
                        regio2.append(prot)
                #print(f'Overlap assigned to Region 2: {regio2}')
                #print(f'------------------- Region 1: {regio1}')
            elif len(regio2) > len(regio1):
                for prot in overlap:
                    if prot in regio2:
                        regio2.remove(prot)
                    if not prot in regio1:
                        regio1.append(prot)
                #print(f'Overlap assigned to Region 1: {regio1}')
                #print(f'------------------- Region 2: {regio2}')
            else:
                decider = random.choice([True,False])
                if decider == True:
                    for prot in overlap:
                        if prot in regio1:
                            regio1.remove(prot)
                        if not prot in regio2:
                            regio2.append(prot)
                    #print(f'Overlap assigned to Region 2: {regio2}')
                    #print(f'------------------- Region 1: {regio1}')
                elif decider == False:
                    for prot in overlap:
                        if prot in regio2:
                            regio2.remove(prot)
                        if not prot in regio1:
                            regio1.append(prot)
                    #print(f'Overlap assigned to Region 1: {regio1}')
                    #print(f'------------------- Region 2: {regio2}')
            
            update_ov_dict = calculate_overlap(dictofcover)
            #print(f'Updated overlaps: {update_ov_dict}')
            #print(f'New cover dict: {dictofcover}')
 
    return dictofcover

    
def summarise_random_walk_results(dictofwalks, nodelist, iterationnum):
    """ 
    Takes in dictionary of walks, nodelist, and iteration number to 
    calculate per protein, of all the visited nodes over all walks,
    how many times each node was visited
    Returns a dictionary with counts
    """
    sampledsubsets = {}
    count = 0
    for prot in nodelist:
        visitednodes = []
        flatlist = []
        for i in range(iterationnum):
            entryname = (prot, count)
            walklist = dictofwalks[entryname]       
            visitednodes.append(walklist)
            count += 1
        for elem in visitednodes:
            if not isinstance(elem, int):
                for subelem in elem:
                    flatlist.append(subelem)
            else:
                flatlist.append(elem)
        count_prots = collections.Counter(flatlist)
        numlist = list(count_prots.most_common())
        #print(numlist)
        #namelist = []
        #for item in numlist:
        #    namelist.append((nodelist[item[0]],item[1]))
        sampledsubsets[prot] = numlist
        #print(f'{prot} {sampledsubsets[prot]}')
    return sampledsubsets

def get_subsets(sampledsubsetdict, nodelist, iteration_num):
    """
    Takes a dictionary of counts and returns subsets per protein 
    in a list, with the constraint that subsets contain proteins 
    visited at least t times 
    """
    t=math.ceil(iteration_num/2)
    masterlist = []
    refdict = {}
    count = 0
    for prot in nodelist:
        refdict[prot] = count
        count += 1

    for prot in nodelist:
        listofsets = []
        listoftups = sampledsubsetdict[prot]
        for pair in listoftups:
            if pair[1] >= t:
                listofsets.append(refdict[pair[0]])
        masterlist.append(set(listofsets))

    return masterlist

def calculate_overlap(dictofsubsets):
    """Takes in a dictionary of subsets, calculates overlaps between all pairwise
    combinations of the subsets, and returns dictionary with set pair IDs and overlaps
    """
    dictofoverlap = {}
    for i in range(len(dictofsubsets)):
        for j in range(i+1, len(dictofsubsets)):
            overlap = list(set(dictofsubsets[i]) & set(dictofsubsets[j]))
            if overlap:
                dictofoverlap[(i,j)] = overlap

    return dictofoverlap   

 
def save_obj(obj, name):
    """Saves object to folder named obj"""
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    """Loads saved object from obj, need to give name"""
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def check_connection(subgraphnodes, originalgraph):
    """ """
    H = originalgraph.subgraph(subgraphnodes)
    val = nx.is_connected(H)

    return val

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s <edgelistfile> <sigfile> <walklength> <iterationnum>",
                                     description="Performs IntCryOmics analysis.")
    parser.add_argument("edgelistfile", help="File with proximity network edges listed", type=str)
    parser.add_argument("significancefile", help="File with significance values per protein", type=str)
    parser.add_argument("walklength", help="Walklength for random walk sampling", type=int)
    parser.add_argument("iterationnum", help="Number of times to repeat random walk", type=int)
    args = parser.parse_args()
    return args

## MAIN
if __name__ == "__main__":
    # parse arguments
    Args = parse_arguments()
    infile = Args.edgelistfile
    sigfile = Args.significancefile
    walk_length = int(Args.walklength)
    iteration_num = int(Args.iterationnum)
    # get respective matrices and print out, get nodes list
    Net, Adjmtx = get_network_from_edges(infile)
    Tmtx = transit_design(Adjmtx)
    Nodelist = list(Net.nodes())

    # iteratively induce random walk, save to a dictionary
    WalksDict = random_walk_new(Net, Tmtx, iteration_num, walk_length)
    # WalksDict = iterate_random_walk(Nodelist, iteration_num, Tmtx, walk_length)

    # read in significance values and check for agreement with nodelist
    #SigDict = get_significance_values(sigfile)
    SigDict = get_sig_vals_full(sigfile)
    Difflist = [x for x in Nodelist if x not in set(list(SigDict.keys()))] 

    # any missing proteins from data are replaced with 0, assumed non-significance
    for missingprot in Difflist:
        SigDict[missingprot] = [0]
  

    # counts visits per node, saves to a dict, keep only proteins that meet threshold cutoff
    SampledSetsDict = summarise_random_walk_results(WalksDict, Nodelist, iteration_num)
    ListofSets = get_subsets(SampledSetsDict, Nodelist, iteration_num)

    # get minimum set cover of the list of sets, save as a dict
    Cover = get_minimum_set_cover(Nodelist, ListofSets)

    print(f'Covering Set: {Cover}\n')
    print('Sampled Regions')
    DictofSets = {}
    ProtNames = []
    for i in range(0, len(Cover)):
        DictofSets[i] = Cover[i] 
        for j in range(0, len(Cover[i])):
            idxlist = list(Cover[i])
            ProtNames.append(Nodelist[idxlist[j]])
        print(f'Region {i}:\t{ProtNames}')
        ProtNames = []
    print('\n')
      

        # calculate overlap between sets in minimum cover
##        DictofOverlap = calculate_overlap(Cover)

        # define regions by randomly reassigning overlaps to shorter regions
##        RegionsDict = remove_overlaps(Cover, DictofOverlap)
##        print('Non-overlapping Regions')
##        for entry in RegionsDict:
##            print(f'{entry} {RegionsDict[entry]} {len(RegionsDict[entry])}')
##            names = []
##            for num in RegionsDict[entry]:
##                names.append(Nodelist[num])
##            Connected = check_connection(names, Net)
##            print(f'{names} {Connected}')
##        print('\n')

    TestDict, ListofPs = fishers_exact_test_regions(DictofSets, Nodelist, SigDict)
    print(TestDict)

    print('\n')
    print(ListofPs)


    #multiple testing correction, bonferroni
    adjust_pvals = multipletests(ListofPs, method='bonferroni')
    print('\nAdjusted p-values after Bonferroni:')
    print(adjust_pvals) 



        # draw the plots according to colormap
#        listofcolors = ["chartreuse","darksalmon","silver","deepskyblue","violet","tomato","limegreen","thistle","lightpink","coral","aliceblue","red","goldenrod"]

#        colormap = ["black" for x in range(len(Nodelist))]
#        count = 0
#        for idxlist in list(RegionsDict.values()):
#            print(f'List of IDs = {idxlist}')
#            for idx in idxlist:
#                colormap[idx] = listofcolors[count] 
#            count += 1

#        print(Net.edges(data=True))
#        elarge = [(u, v) for (u, v, d) in Net.edges(data=True) if d['weight'] > 0.5]
#        esmall = [(u, v) for (u, v, d) in Net.edges(data=True) if d['weight'] <= 0.5]
#        pos = nx.spring_layout(Net)
#        nx.draw_networkx_nodes(Net, pos, node_size=500, node_color=colormap)
#        nx.draw_networkx_edges(Net, pos, edgelist=elarge, width=6)
#        nx.draw_networkx_edges(Net, pos, edgelist=esmall, width=6, alpha=0.5, edge_color='b', style='dashed')
#        nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')
##        nx.draw(Net,node_size=500, node_color=colormap,with_labels=True)
    plt.show()

        # test regions, binomial and hypergeometric
        #BinomDict = binom_test_regions(RegionsDict, Nodelist, SigDict, siglvl=0.05)
        #print(BinomDict)
        #for entry in BinomDict:
        #    print(f'{entry} {BinomDict[entry]}')


