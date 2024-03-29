#!/usr/bin/env python3
"""
[intcryomics_sigassign.py]

This script is a command-line interface to perform the main COSNet_i analysis.

It creates a structural network, samples regions, and tests each region
for significance according to given omics data.

You have the option to assign your own significance values, or upload a significance file.

USAGE: python3 intcryomics.py filewithedgelist walklength iterationnum --sigfile [sigfile]
"""
from sys import exit
import collections
import os
import time
import random
import math
import pickle
import itertools
import argparse
from pathlib import Path

from scipy import stats
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests

def get_network_from_edges(edgelistfile):
    """ 
    Creates network from list of connected nodes.
    Gets adjacency matrix.

    Parameter
    ---------
    edgelistfile: pathlib.PosixPath

    Returns
    -------
    G: networkx.classes.graph.Graph
    adjmtx: numpy.ndarray
    """
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

    Parameter
    ---------
    adjacent_mtx: numpy.ndarray

    Returns
    -------
    transmtx: numpy.ndarray
    """
    summtx = np.diag(np.sum(adjacent_mtx, axis=0))
    transmtx = np.dot(np.linalg.inv(summtx),adjacent_mtx)
    np.savetxt('Tmatrix.dat', transmtx, delimiter = ',')    
    return transmtx

def generateSequence(startIndex, transitionMatrix, path_length, alpha = 0.1):
    """ 
    Generates random sequence of steps according to the transition matrix.

    Parameters
    ----------
    startIndex: int
        index of starting node
    transitionMatrix: numpy.ndarray
        matrix of transition probabilities
    path_length: int
    alpha: float
        allows for circular steps

    Returns
    -------
    result: list
        list of steps in random sequence
    """
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
    """
    Performs random walk for given length with number
    of repeats given by iter_num.

    Parameters
    ----------
    G: networkx.classes.graph.Graph
    transitmtx: numpy.ndarray
    iter_num: int
    walk_length: int

    Returns
    -------
    WalksDict: dict
        dict of walks of format {(startnode, walk#): [step1, step2, ...], ...}
    """
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

    Parameters
    ----------
    filewithsigs: pathlib.PosixPath
    
    Returns
    -------
    sig_dict: dict
        dict of form: { prot : 0/1 }
    """
    sig_dict = {}
    with open(filewithsigs, 'r') as f:
        for line in f:
            line_elems = line.strip().split(' ')
            sig_dict[line_elems[0]] = int(line_elems[1])
    return sig_dict

def manual_get_sig_vals(regiondict, listofprots):
    """
    Allows for manual assignment of significance values.

    Parameters
    ----------
    regiondict: dict
    listofprots: list

    Returns
    -------
    altsigdict: dict
    """
    altsigdict={}
    for prot in listofprots: #set all prots to 0
        altsigdict[prot]=[0] 

    proceed=input(f"\nYou have not given a significance file.\n\
You will have to enter significance values per protein individually. Default sig values are 0.\n\
Will you manually enter significance values? [y/n] ")
    if proceed in 'yYesyesYES':
        regionstoassign=input("Refer to the regions defined above.\n\
Type in the region number(s) in which you would like to assign significances (separated by one space): ")
        regs=[int(item) for item in regionstoassign.split()]
        for region_nr in regs:
            print('\n')
            randomassign=input(f'---> Assigning sig values to Region {region_nr}.\n\
     You can do this randomly or manually to each protein. Do you want to do it randomly? [y/n] ')
            if randomassign in 'y Yes yes YES':
                numprots = len(regiondict[region_nr])
                print(f'     You have {numprots} proteins in Region {region_nr}.')
                numsig=int(input(f'     Enter # significances to randomly assign: '))
                if int(numsig) <= numprots:
                    randlist=random.sample(range(0,numprots),numsig)
                    if randlist:
                        for idx in randlist:
                            protname=regiondict[region_nr][idx]
                            print(f'         -Assigning 1 to {protname}')
                            altsigdict[protname]=[1]
                            time.sleep(0.2)
                else:
                    print('Please enter a valid number. Start again...')
                    break
            else:
                print(f'     You have chosen to assign significances manually.\n     Here are the proteins in Region {region_nr}:')
                for i in range(0,len(regiondict[region_nr])):
                    print(f'     {i}: {regiondict[region_nr][i]}')
                sigidx=input("     Please enter the indices corresponding to proteins you want to assign significant (space-separated):\n   ->") 
                sigidxlst=[int(item) for item in sigidx.split()]
                for protidx in sigidxlst:
                    protname=regiondict[region_nr][protidx]
                    print(f'         -Assigning 1 to {protname}')
                    altsigdict[protname]=[1]
                    time.sleep(0.2)
                    
        print('\nManual assignment of significances completed. Continuing with analysis...\n')
        time.sleep(0.6)
    else:
        altfile=input('Enter the path to and name of a sig file, or Q to quit: ')
        if altfile in 'Q qQuitquit':
            exit('User terminated. Please rerun if necessary.')
        elif Path(altfile).is_file():
            altsigdict=get_sig_vals_full(Path(altfile))

    return altsigdict

def get_sig_vals_full(filewithsigs):
    """
    Gets dict of omics significance values per node in network.
    Gives back form of {uL1:[0,0,0],uL2:[1,1,1]...}

    Parameters
    ----------
    filewithsigs: pathlib.PosixPath
        whitespace separated, two col file with nodes and significances

    Returns
    -------
    sig_dict_full: dict
        dict of nodes and a list of their significances
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
    Performs fishers exact test for significance per sampled region.

    Parameters
    ----------
    dictofregions: dict
    Nodelist: list
    SigDict: dict

    Returns
    -------
    TestResDict: dict
    PvalList: list
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

    Parameters
    ----------
    nodelist: list
    listofsubsets: list

    Returns
    -------
    cover: list
        list of sets of sampled regions
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
    how many times each node was visited.
    Returns a dictionary with counts.

    Parameters
    ----------
    dictofwalks: dict
    nodelist: list
    iterationnum: int

    Returns
    -------
    sampledsubsets: dict
        dict of nodes visited and visitation counts per node
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
    visited at least t times.

    Parameters
    ----------
    sampledsubsetdict: dict
        dict of nodes visited and visitation counts per node
    nodelist: list
    iteration_num: int

    Returns
    -------
    masterlist: list
        subsets per protein
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

    Parameters
    ----------
    dictofsubsets: dict

    Returns
    -------
    dictofoverlap: dict
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

def write_out_regions(namesdict):
    regionsfile='sampled_regions.txt'
    with open(regionsfile, 'w+') as f:
        for k, v in namesdict.items():
            f.write(f'Region {k}: {v}\n')
    return

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s edgelistfile walklength iterationnum --sigfile [significancefile]",
                                     description="Performs IntCryOmics analysis.")
    parser.add_argument("edgelistfile", help="File with proximity network edges listed", type=str)
    parser.add_argument("walklength", help="Walklength for random walk sampling", type=int)
    parser.add_argument("iterationnum", help="Number of times to repeat random walk", type=int)
    parser.add_argument("--sigfile", help="Optional file with significance values per protein", nargs=1, type=str)
    args = parser.parse_args()
    return args

## MAIN
if __name__ == "__main__":
    # parse arguments
    Args = parse_arguments()
    infile = Args.edgelistfile
    walk_length = Args.walklength
    iteration_num = Args.iterationnum
    # get respective matrices and print out, get nodes list
    Net, Adjmtx = get_network_from_edges(infile)
    Tmtx = transit_design(Adjmtx)
    Nodelist = list(Net.nodes())


    # iteratively induce random walk, save to a dictionary
    WalksDict = random_walk_new(Net, Tmtx, iteration_num, walk_length)

    # counts visits per node, saves to a dict, keep only proteins that meet threshold cutoff
    SampledSetsDict = summarise_random_walk_results(WalksDict, Nodelist, iteration_num)
    ListofSets = get_subsets(SampledSetsDict, Nodelist, iteration_num)

    # get minimum set cover of the list of sets, save as a dict
    Cover = get_minimum_set_cover(Nodelist, ListofSets)

    print(f'Covering Set: {Cover}\n')
    print('Sampled Regions')
    DictofSets = {}
    DictofNames = {}
    ProtNames = []
    for i in range(0, len(Cover)):
        DictofSets[i] = Cover[i] 
        for j in range(0, len(Cover[i])):
            idxlist = list(Cover[i])
            ProtNames.append(Nodelist[idxlist[j]])
        DictofNames[i]=ProtNames
        print(f'Region {i}:\t{ProtNames}')
        ProtNames = []
    print('\n')
  
    # write out regions to a file
    write_out_regions(DictofNames)
    print('\n---> Regions written to sampled_regions.txt')

    # read in significance values and check for agreement with nodelist
    #SigDict = get_significance_values(sigfile)
    if Args.sigfile:
        SigDict = get_sig_vals_full(Path(Args.sigfile[0]))
        Difflist = [x for x in Nodelist if x not in set(list(SigDict.keys()))] 

        # any missing proteins from data are replaced with 0, assumed non-significance
        for missingprot in Difflist:
            SigDict[missingprot] = [0]
    else:
        SigDict=manual_get_sig_vals(DictofNames, Nodelist)
        if type(SigDict) is dict and SigDict: 
            Difflist = [x for x in Nodelist if x not in set(list(SigDict.keys()))] 
            for missingprot in Difflist:
               SigDict[missingprot] = [0]
        else:
            print('WARNING: Manual assignment unsuccessful.')
            exit('Terminated. Please rerun.')
  
    # fishers exact test 
    TestDict, ListofPs = fishers_exact_test_regions(DictofSets, Nodelist, SigDict)
    print(TestDict)

    print('\n')
    print(ListofPs)


    #multiple testing correction, bonferroni
    adjust_pvals = multipletests(ListofPs, method='bonferroni')
    print('\nAdjusted p-values after Bonferroni:')
    print(adjust_pvals) 

