#!usr/env/bin/ python3
"""
combination.py: takes in a list of file names, 
makes another list with all dual combinations

python3 combination.py [listoffiles] [combnum] [outputfile] [outpath]

"""
## IMPORTS
from sys import argv
import itertools
import os
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
    """
    """
    filename = os.path.join(outpath,filename)
    outfile = open(filename, 'w+')
    for line in combinedlist:
        ent1 = line[0]
        ent2 = line[1]
        outfile.write(f'{ent1} {ent2}\n')
    outfile.close()

## MAIN
if __name__ == "__main__":
    if len(argv) != 5:
        print('################################################################################################\n')
        print('    Usage: python3 combination.py [filewithlistoffiles] [combi_num] [outfilename] [outpath]     \n')
        print('################################################################################################\n')
    else:
        inputfile = argv[1]
        combnum = int(argv[2])
        outputfile = argv[3]
        outpath = argv[4]
        FileList = get_list(inputfile)
        CombList = combinate(FileList,combnum)
        print_output(CombList,outputfile,outpath)
