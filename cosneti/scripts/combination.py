#!usr/env/bin/ python3
"""
combination.py: 

Reads in file with list of file names,  
Outputs another file  with all n-combinations of original files,
one combination per line.

Useful for preparing protein pdb file input for pair-wise contact calculation.

python3 combination.py [listoffiles] [combnum] [outputfile] [outpath]

"""
## IMPORTS
from sys import argv
import itertools
import os
import argparse
## FUNCTIONS
def get_list(infile):
    """
    Puts lines in file in a list
    """
    filelist = []
    with open(infile) as f:
        for line in f:
            filelist.append(line.strip())

    return filelist

def combinate(listoffiles, combnum):
    """
    Makes a list of all the N choose combnum options
    """
    combinedlist = list(itertools.combinations(listoffiles, combnum))
    return combinedlist

def print_output(combinedlist, filename, outpath):
    filename = os.path.join(outpath,filename)
    outfile = open(filename, 'w+')
    for item in combinedlist:
        ent = " ".join(item)
        outfile.write(f'{ent}\n')
    outfile.close()

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] <inputfile> <N> <outfile> <outpath",
                                     description="Returns file of N-wise combinations of input file list.")
    parser.add_argument("inputfile", help="One column text file of filenames", type=str)
    parser.add_argument("N", help="Combination number, e.g. 2 gives pair-wise combinations", type=int)
    parser.add_argument("outfile", help="Name of outputfile", type=str)
    parser.add_argument("outpath", help="Path to output outfile", type=str)
    args = parser.parse_args()
    return args

def main():
    Args = parse_arguments()
    FileList = get_list(Args.inputfile)
    if Args.N == 1:
        print(f'Combination of 1 is just your input file: {Args.inputfile}')
    elif Args.N > len(FileList):
        print(f'Combination number: {Args.N} cannot be larger than number of files inputted: {len(FileList)}')
    else:
        CombList = combinate(FileList,Args.N)
        print_output(CombList,Args.outfile,Args.outpath)
    
## MAIN
if __name__=="__main__":
    main()
