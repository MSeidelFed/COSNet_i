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
import numpy as np
import os
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
    ent_id = struct_id_elems[1].split(".")[0]

    structure = parser.get_structure(ent_id, PDBfilepath)
    residues = structure.get_residues()
    atoms = structure.get_atoms()

    return structure, struct_id, ent_id, residues

def center_of_mass(one_residue, flag):
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

## MAIN
if __name__ == "__main__":
    #usage
    if len(argv) != 4:
        print('################################################################################################\n')
        print('     Usage: python3 calculate_distance.py [pathtoPDBfile1] [pathtoPDBfile2] [outputpath]        \n')
        print('################################################################################################\n')
    else:
        #read in comparison files
        ent1_file = argv[1]
        ent2_file = argv[2]
        outputpath = argv[3]
        if outputpath == ".":
            outputpath = os.getcwd()
        coms_method = 0

        #parse out structures and info of both files
        Ent1Struct, Struct1ID, Ent1ID, Residues1 = get_struct(ent1_file)
        Ent2Struct, Struct2ID, Ent2ID, Residues2 = get_struct(ent2_file)

        #calculate center of masses for all residues
        ComsList1 = []
        for res in Residues1:
            com1 = center_of_mass(res, coms_method)
            ComsList1.append(com1)

        ComsList2 = []
        for res in Residues2:
            com2 = center_of_mass(res, coms_method)
            ComsList2.append(com2)

        DistMtx = calc_dist_matrix(ComsList1, ComsList2)
    
        #save distance matrix into output file
        if Struct1ID != Struct2ID:
            np.savetxt('dist_mtx_{}_{}.csv'.format(Struct1ID,Struct2ID), DistMtx, delimiter=",")
            outfile = f'dist_mtx_{Struct1ID}_{Struct2ID}.csv'
            os.rename(outfile, os.path.join(outputpath, outfile))
        else:
            np.savetxt('dist_mtx_{}_{}.csv'.format(Ent1ID,Ent2ID), DistMtx, delimiter=",")
            outfile = f'dist_mtx_{Ent1ID}_{Ent2ID}.csv'
            os.rename(outfile, os.path.join(outputpath, outfile))
            
