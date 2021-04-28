#!/usr/bin/env python3
"""
[batch_calc_dist.py]

This script is a command-line interface to calculate distances
between structures in batch. It runs calculate_distance.py over
a given list of pdb files pairs.

USAGE: python3 batch_calc_dist.py infile inputpath outputpath
"""
import argparse
import subprocess
from pathlib import Path

def iterate_run_calcdist(listoffiletuples, pdbpath, outpath):
    """
    Iteratively runs calculate_distance.py over pairs of files.

    Parameters
    ----------
    listofffiletuples: list
    pdbpath: pathlib.PosixPath
    outpath: pathlib.PosixPath
    """
    counter=1
    for filetuple in listoffiletuples:
        pdbfile1 = pdbpath.joinpath(filetuple[0])
        pdbfile2 = pdbpath.joinpath(filetuple[1])
        print(f"Running on protein # {counter}")
        print(f"{pdbfile1}")
        print(f"{pdbfile2}")
        subprocess.run(["python3", "calculate_distance.py", f"{pdbfile1}", f"{pdbfile2}", f"{outpath}"])
        counter+=1

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] infile inputpath outputpath",
                                     description="Runs calculate_distance.py iteratively over combinations of pdb files")
    parser.add_argument("infile", help="file with dual combinations of pdb files from combination.py")
    parser.add_argument("inputpath", help="path to where pdb files are found")
    parser.add_argument("outputpath", help="path to output distance matrices, needed for calculate_distance.py")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    Args = parse_arguments()
    infile = Path(Args.infile)
    inputpath = Path(Args.inputpath)
    outputpath = Path(Args.outputpath)
    if infile.is_file():
        runlist=[]
        with open(infile, 'r') as f:
            for line in f.readlines():
                runlist.append((line.strip().split()[0], line.strip().split()[1])) 
        if runlist:
            iterate_run_calcdist(runlist, inputpath, outputpath)
    else:
        raise FileNotFoundError('Give a correct combinated file!')
