#!/usr/env/bin/ python3
"""
[split_cif_by_entity.py]

USAGE: python3 split_cif_by_entity.py <inputfile> <pathtooutput>
"""
###### IMPORTS ######
from sys import argv
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from collections import Counter
import argparse
import os
import glob
## FUNCTIONS
def lower_upper_idx(EntityAtomIDs):
    """
    Returns dictionary of lower and upper indices for 
    entities in the atomic data block of mmCIF file,
    """
    ent_instances = dict(Counter(EntityAtomIDs))
    idx_dict = {}
    for i in range(len(ent_instances)):
        if i == 0:
            end = ent_instances['1']-1
            idx_dict[i] = (0,end)
        elif i != 0:
            length = ent_instances[str(i+1)]
            idx_dict[i] = (end+1,end+length)
            end = end+length
    return idx_dict

def grab_struct_data(ObjName, IDxDict, NamesList, CIFDict, address):
    """
    Takes values from cifdict corresponding to the desired
    structural features, prints them according to PDB standard
    format in a file (column ranges serve as identifiers)
    """
    if address[:-1] != "/":
        address = address + "/"

    atomlist = CIFDict["_atom_site.group_PDB"]
    atomid_num = CIFDict["_atom_site.id"]
    atomid_lbl = CIFDict["_atom_site.label_atom_id"]
    atomtype = CIFDict["_atom_site.type_symbol"]
    atom_comp = CIFDict["_atom_site.label_comp_id"]
    atom_asym = CIFDict["_atom_site.label_asym_id"]
    atom_seqid = CIFDict["_atom_site.label_seq_id"]
    
    atom_Xcoord = CIFDict["_atom_site.Cartn_x"]
    atom_Ycoord = CIFDict["_atom_site.Cartn_y"]
    atom_Zcoord = CIFDict["_atom_site.Cartn_z"]

    atom_occupancy = CIFDict["_atom_site.occupancy"]
    atom_iso = CIFDict["_atom_site.B_iso_or_equiv"]

    for i in range(len(IDxDict)):
        entity_name = NamesList[i].replace(" ","_")
        print(f'{i}: Creating pdb file for {entity_name}\n')
        if 'protein' in entity_name: 
            protid = entity_name.split('_protein_')
            entity_name = f'{protid[1]}_{protid[0]}'
        structfile = open(f'{address}{ObjName}_{entity_name}.pdb', 'w+')
        bounds = IDxDict[i]
        k = 1
        for j in range(bounds[0],bounds[1]+1):
            structfile.write('{0:<6}{1:>5} {2:<4} {3:<3}{4:>2}{5:>4}    {6:>8}{7:>8}{8:>8}{9:>6}{10:>6}          {11:>2}'.format(
                            atomlist[j],k,atomid_lbl[j],atom_comp[j],atom_asym[j],
                            atom_seqid[j],atom_Xcoord[j],atom_Ycoord[j],atom_Zcoord[j],
                            atom_occupancy[j],atom_iso[j],atomtype[j]))
            k += 1
            if j != bounds[1]:
                structfile.write('\n')
        structfile.close()

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] <pathto/CIFfile> <outpath>",
                                     description="Splits CIF into separate PDB files per entity")
    parser.add_argument("CIFfile", help="CIF file input to be split", type=str)
    parser.add_argument("outpath",help="Output path", type=str)
    args = parser.parse_args()
    return args
    
## MAIN
if __name__ == "__main__":
    # get arguments
    Args = parse_arguments()
    # read in data file from command line
    pathtociffile = Args.CIFfile
    obj_name = os.path.basename(pathtociffile)[:-4]    
    address = Args.outpath
    # save contents as a dictionary
    cif_dict = MMCIF2Dict(pathtociffile)
    # get lists of identifiers and entities
    ent_id = cif_dict["_entity_poly.entity_id"]
    ent_name = cif_dict["_entity.pdbx_description"]
    ent_atom_ids = cif_dict["_atom_site.label_entity_id"]
    # get dictionary of indices per entity indicating corresponding rows
    # of the atomic data 
    IDxDict = lower_upper_idx(ent_atom_ids)
    grab_struct_data(obj_name, IDxDict, ent_name, cif_dict, address)

