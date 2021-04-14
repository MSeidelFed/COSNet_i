#!/usr/env/bin python3
"""
[batch_reindex_pdb.py]

DESCRIPTION

USAGE

RETURNS

"""
### IMPORTS
import argparse
import subprocess
from pathlib import Path

### FUNCTIONS
def iterate_run_reindex(listofpdbtuples, pdbpath, outpath):
    counter=1
    for filetuple in listofpdbtuples:
        pdbfile = pdbpath.joinpath(filetuple[0])
        startidx = filetuple[1]
        outfile = f"{outpath.joinpath(filetuple[2])}"
        print(f"Running on protein # {counter}")
        print(f"Input: {pdbfile} Output: {outfile}")
        print(f"Start index {startidx}")
        subprocess.run(["python3", "reindex_pdb.py", f"{startidx}", f"{pdbfile}", f"{outfile}"])
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
