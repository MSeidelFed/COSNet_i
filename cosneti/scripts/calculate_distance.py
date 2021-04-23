#!/usr/env/bin/ python3
"""
[calculate_distance.py]

DESCRIPTION: Calculates distances between objects in two entities 
(e.g. residues in PDB files) and stores in a distance matrix

USAGE: python3 calculate-distance.py <pdbfile1> <pdbfile2>

RETURNS: dist_mtx_<pdbfile1>_<pdbfile1>.csv

Currently uses the geometric center of mass for coarse-graining 
the entity (a protein amino acid represented by the center of mass
of its atoms, for example, and same goes for a nucleotide)
"""
###### IMPORTS ######
from sys import argv
from Bio.PDB import PDBParser
from pathlib import Path
import numpy as np
import os
import argparse
## FUNCTIONS
def get_struct(PDBfilepath):
    """
    Takes in PDB file, parses out structural information
    Returns structure object, original PDB ID of structure, entity name (if 
    file is a portion of the original larger file, say an r protein or rRNA) 
    """
    parser = PDBParser(PERMISSIVE=1)
    PDBfile = os.path.basename(PDBfilepath)
    struct_id_elems = PDBfile.split("_")
    struct_id = struct_id_elems[0]
    if len(struct_id_elems)!= 2:
        ent_id = PDBfile.split('.')[0]
    else:
        ent_id = struct_id_elems[1].split(".")[0]

    structure = parser.get_structure(ent_id, PDBfilepath)
    residues = structure.get_residues()
    atoms = structure.get_atoms()

    return structure, struct_id, ent_id, residues

def center_of_mass(one_residue, flag=0):
    """
    Finds the center of mass of a residue
    Returns this as (x,y,z)
    """
    x, y, z = 0, 0, 0
    com = 0
    if flag == 0:
        for atom in one_residue:
            res_size = len(one_residue)
            x += atom.coord[0]
            y += atom.coord[1]
            z += atom.coord[2]
        com = (x/res_size, y/res_size, z/res_size)
        return com
    else:
    # here do the gravitational-based centre of mass  
        return com
      
def calc_resi_dist(coords1, coords2):
    """
    Takes two coordinates, calculates euclidean distance
    between the two
    """
    a = np.array(coords1)
    b = np.array(coords2)
    dist = np.linalg.norm(a-b)
    return dist

def calc_dist_matrix(entity1, entity2):
    """
    Takes in two lists of coordinates, creates distance
    matrix to store all pairwise combinations of the coordinates
    """
    dist_mtx = np.zeros((len(entity1),len(entity2)))
    for i, coords_one in enumerate(entity1):
        for j, coords_two in enumerate(entity2):
            dist_mtx[i,j] = calc_resi_dist(coords_one, coords_two)
    return dist_mtx

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] <pathto/PDBfile1> <pathto/PDBfile2> <outpath>",
                                     description="Calculates Euclidean distances between residues of two PDB files, stores in distance matrix file.")
    parser.add_argument("pdbfile1",help="First PDB file",type=str)
    parser.add_argument("pdbfile2",help="Second PDB file",type=str)
    parser.add_argument("outpath",help="Where to output distance matrix file",type=str)
    args = parser.parse_args()
    return args

## MAIN
if __name__ == "__main__":
    #usage
    Args = parse_arguments()
    if Args.outpath=='.':
        outpath=Path.cwd()
    else:
        outpath=Path(Args.outpath)
    #parse out structures and info of both files
    Ent1Struct, Struct1ID, Ent1ID, Residues1 = get_struct(Args.pdbfile1)
    Ent2Struct, Struct2ID, Ent2ID, Residues2 = get_struct(Args.pdbfile2)

    #calculate center of masses for all residues
    ComsList1 = []
    for res in Residues1:
        com1 = center_of_mass(res)
        ComsList1.append(com1)

    ComsList2 = []
    for res in Residues2:
        com2 = center_of_mass(res)
        ComsList2.append(com2)

    DistMtx = calc_dist_matrix(ComsList1, ComsList2)
    
    #save distance matrix into output file
    if Struct1ID != Struct2ID:
        np.savetxt('dist_mtx_{}_{}.csv'.format(Struct1ID,Struct2ID), DistMtx, delimiter=",")
        outfile = f'dist_mtx_{Struct1ID}_{Struct2ID}.csv'
        os.rename(outfile, os.path.join(outpath, outfile))
    else:
        np.savetxt('dist_mtx_{}_{}.csv'.format(Ent1ID,Ent2ID), DistMtx, delimiter=",")
        outfile = f'dist_mtx_{Ent1ID}_{Ent2ID}.csv'
        os.rename(outfile, os.path.join(outpath, outfile))