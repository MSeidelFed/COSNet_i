#!/usr/env/bin python3
"""
[batch_calc_dist.py]

DESCRIPTION

USAGE

RETURNS

"""
### IMPORTS
import glob
import sys
import os
import argparse
import subprocess
from pathlib import Path

### PATH TO SCRIPTS
sys.path.append(".")
sys.path.append("Python_Modules")
sys.path.append("../Python_Modules")

### FUNCTIONS
def iterate_run_calcdist(listoffiletuples, pdbpath, outpath):
    counter=1
    for filetuple in listoffiletuples:
        pdbfile1 = pdbpath.joinpath(filetuple[0])
        pdbfile2 = pdbpath.joinpath(filetuple[1])
        print(f"Running on protein # {counter}")
        print(f"{pdbfile1}")
        print(f"{pdbfile2}")
        subprocess.run(["python3", "calculate_distance.py", f"{pdbfile1}", f"{pdbfile2}", f"{outpath}"])
        counter+=1
    return

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] infile inputpath outputpath",
                                     description="Runs calculate_distance.py iteratively over combinations of pdb files")
    parser.add_argument("infile", help="file with dual combinations of pdb files from combination.py")
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
                runlist.append((line.strip().split()[0], line.strip().split()[1])) 
        if runlist:
            iterate_run_calcdist(runlist, inputpath, outputpath)
    else:
        raise FileNotFoundError('Give a correct combinated file!')
