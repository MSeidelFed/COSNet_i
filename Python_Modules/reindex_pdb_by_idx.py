#!/usr/bin/env python3
"""
[reindex_pdb_residues.py]

Returns a reindexed pdb file from a given starting index.
"""
### IMPORTS
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from pathlib import Path
import argparse

### FUNCTIONS
def parse_pdb(PDBfile):
    """Parses PDB file, returns PDB object"""
    name=PDBfile.name[:-4]
    p=PDBParser()
    struct = p.get_structure(name, PDBfile)
    return struct, name

def reindex(structobj, startidx):
    """Reindexes residues per chain per model"""
    for model in structobj: 
        for chain in model:
            count = startidx
            for resi in chain.get_residues():
                if resi.id[0] == ' ': 
                    resi._id = (' ', count, ' ')
                    count += 1
    return structobj

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] startidx pdbfile outfile",
                                     description="Returns PDB file with residues reindexed according to start index")
    parser.add_argument("startidx", help="Starting index, usually 1", type=int)
    parser.add_argument("pdbfile", help="PDB file to be reindexed", type=str)
    parser.add_argument("outfile", help="Reindexed PDB file", type=str)
    args = parser.parse_args()
    return args
    
def main():
    # Get arguments
    Args = parse_arguments()
    # Reindex PDB structure
    p = Path(Args.pdbfile)
    o = Path(Args.outfile)
    if p.is_file() and o.parent.is_dir():
        Struct, Name = parse_pdb(p)
        Reindexed = reindex(Struct, Args.startidx)
        # Save reindexed structure
        io = PDBIO()
        io.set_structure(Reindexed)
        io.save(Args.outfile)
        print(f'Reindexed file: {Name}_reindex.pdb')
        print(f'Saved at this location: {o.parent}')
    else:
        print("Problem with input file or output path, check that both exist!")

### MAIN
if __name__=="__main__":
    main()
