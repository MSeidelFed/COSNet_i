#!/usr/bin/env python3
"""
[combination.py] 

USAGE: python3 combination.py listoffiles combnum outputfile outpath

This script is a command-line interface to make n-combinations of 
files. It is useful for preparing protein pdb file input for pair-wise
contact calculation. 
"""

from sys import argv
from pathlib import Path
import itertools
import argparse

def get_list_of_files(infile):
    """
    Makes list of file names 

    Parameters
    ----------
    infile: pathlib.PosixPath

    Returns
    -------
    filelist: list
        list of file names
    """
    filelist = []
    with open(infile) as f:
        for line in f:
            filelist.append(line.strip())

    return filelist

def combinate(listoffiles, combnum):
    """
    Makes a list of tuples with the N choose combnum options

    Parameters
    ----------
    listoffiles: list
        list of file names
    combnum: int
        combination number (N choose combnum)

    Returns
    -------
    combinedlist: list
        list of tuples with combnum-combinations of files
    """
    combinedlist = list(itertools.combinations(listoffiles, combnum))
    return combinedlist

def print_output(combinedlist, filename, outpath):
    """
    Writes out file with combinations of files per line

    Parameters
    ----------
    combinedlist: list
        list of tuples with file combinations
    filename: str
        name of file to write file combinations
    outpath: pathlib.PosixPath
        path to output file
    """
    filename = outpath.joinpath(filename)
    outfile = open(filename, 'w+')
    for item in combinedlist:
        ent = " ".join(item)
        outfile.write(f'{ent}\n')
    outfile.close()

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] inputfile N outfile outpath",
                                     description="Returns file of N-wise combinations of input file list.")
    parser.add_argument("inputfile", help="One column text file of filenames", type=str)
    parser.add_argument("N", help="Combination number, e.g. 2 gives pair-wise combinations", type=int)
    parser.add_argument("outfile", help="Name of outputfile", type=str)
    parser.add_argument("outpath", help="Path to output outfile", type=str)
    args = parser.parse_args()
    return args

def main():
    Args = parse_arguments()
    inputfile = Path(Args.inputfile)
    if not inputfile.is_file():
        raise FileNotFoundError('Check your input file!')
    else:
        FileList = get_list_of_files(inputfile)
        if Args.N == 1:
            print(f'Combination of 1 is just your input file: {Args.inputfile}')
        elif Args.N > len(FileList):
            print(f'Combination number: {Args.N} cannot be larger than number of files inputted: {len(FileList)}')
        else:
            CombList = combinate(FileList,Args.N)
            print_output(CombList,Path(Args.outfile),Path(Args.outpath))
    
if __name__=="__main__":
    main()
