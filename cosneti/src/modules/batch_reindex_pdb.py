#!/usr/env/bin python3
"""
[batch_reindex_pdb.py]

This script is a command-line interface to reindex pdb files in batch.
It runs reindex_pdb.py over a given list of pdbfiles.

USAGE: python3 batch_reindex_pdb.py infile inputpath outputpath
"""
### IMPORTS
import sys
import argparse
from pathlib import Path

### PATH TO SCRIPTS
sys.path.append(".")
sys.path.append("Python_Modules")
sys.path.append("../Python_Modules")

from reindex_pdb import *

### FUNCTIONS
def iterate_run_reindex(listofpdbtuples, pdbpath, outpath):
    """
    Iteratively runs reindex_pdb.py over pdb files.

    Parameters
    ----------
    listofpdbtuples: list
    pdbpath: pathlib.PosixPath
    outpath: pathlib.PosixPath
    """
    counter=1
    for filetuple in listofpdbtuples:
        pdbfile = pdbpath.joinpath(filetuple[0])
        startidx = filetuple[1]
        outfile = f"{outpath.joinpath(filetuple[2])}"
        print(f"Running on protein # {counter}")
        print(f"Input: {pdbfile} Output: {outfile}")
        print(f"Start index {startidx}")
        reindex_pdb(int(startidx), pdbfile, Path(outfile)) 
        counter+=1
    return

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] infile inputpath outputpath",
                                     description="Runs reindex_pdb.py iteratively over given list of pdb files")
    parser.add_argument("infile", help="3-column file with pdb files, start idx, and outputfile names")
    parser.add_argument("inputpath", help="path to where pdb files are found")
    parser.add_argument("outputpath", help="path to output distance matrices, needed for calculate_distance.py")
    args = parser.parse_args()
    return args

### MAIN
if __name__=="__main__":
    Args = parse_arguments()
    infile = Path(Args.infile)
    inputpath = Path(Args.inputpath)
    outputpath = Path(Args.outputpath)
    if infile.is_file():
        runlist=[]
        with open(infile, 'r') as f:
            for line in f.readlines():
                contents=line.strip().split()
                if not len(contents) == 3:
                    runlist.append((contents[0], 1, contents[0][:-4]+'_reindex.pdb'))
                else:
                    runlist.append((contents[0], int(contents[1]), contents[2]))
        if runlist:
            iterate_run_reindex(runlist, inputpath, outputpath)
    else:
        raise FileNotFoundError('Give a correct pdb file!')
